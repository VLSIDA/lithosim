#!/usr/bin/env python3
"""
Lithography simulator with GDS support and physically-based models.

Optical model: Hopkins partially-coherent imaging via SOCS decomposition.
Resist model:  Lumped parameter (Gaussian diffusion blur + threshold).
Input:         GDS/OASIS (via gdstk) or bitmap images (PBM/PNG).
Output:        Aerial images, resist images, binary contours.

Usage:
    python lithosim.py [options] <input> <output_prefix>

Examples:
    python lithosim.py layout.gds out/sim
    python lithosim.py --layer 1/0 --nm-per-pixel 4 chip.gds results/layer1
    python lithosim.py --sigma 0.7 mask.pbm results/test
    python lithosim.py --source annular --sigma-inner 0.5 --sigma-outer 0.9 layout.gds out
"""

import argparse
import sys
import time
import os

import numpy as np
from scipy.special import j1 as bessj1
from scipy.signal import fftconvolve
from scipy.ndimage import gaussian_filter
from skimage import measure
from PIL import Image


# ======================================================================
# GDS / mask loading
# ======================================================================

def load_gds(filename, layer=None, nm_per_pixel=1.0, padding=0):
    """Rasterize a GDS/OASIS layout to a binary mask array.

    Parameters
    ----------
    filename : str
        Path to .gds or .oas file.
    layer : tuple (layer, datatype) or None
        Which layer to rasterize.  If None, uses the first layer found.
    nm_per_pixel : float
        Grid resolution in nm.  GDS coordinates are in um, so
        1 GDS unit = 1000 nm by default (assumes 1um database unit).
    padding : int
        Extra pixels of padding around the bounding box.

    Returns
    -------
    mask : ndarray (float64, 2D)
        1.0 where there is geometry, 0.0 elsewhere.
    origin : tuple (x_nm, y_nm)
        World coordinates of pixel (0,0) for back-annotation.
    """
    import gdstk

    lib = gdstk.read_gds(filename) if filename.endswith('.gds') else gdstk.read_oas(filename)
    top_cells = lib.top_level()
    if not top_cells:
        raise ValueError(f"No top-level cells in {filename}")
    cell = top_cells[0]

    # collect polygons, optionally filtering by layer
    all_polys = cell.get_polygons()
    if layer is not None:
        polys = [p for p in all_polys if p.layer == layer[0] and p.datatype == layer[1]]
    else:
        polys = all_polys

    if not polys:
        layers_present = sorted(set((p.layer, p.datatype) for p in all_polys))
        raise ValueError(f"No polygons found on layer {layer}. Available: {layers_present}")

    # database unit: GDS stores in um, convert to nm
    db_unit_nm = lib.unit * 1e9  # lib.unit is in meters

    # compute bounding box in nm
    all_points = np.vstack([p.points for p in polys]) * db_unit_nm
    x_min, y_min = all_points.min(axis=0) - padding * nm_per_pixel
    x_max, y_max = all_points.max(axis=0) + padding * nm_per_pixel

    width = int(np.ceil((x_max - x_min) / nm_per_pixel))
    height = int(np.ceil((y_max - y_min) / nm_per_pixel))
    mask = np.zeros((height, width), dtype=np.float64)

    # rasterize each polygon using scanline fill
    for poly in polys:
        pts_nm = poly.points * db_unit_nm
        # convert to pixel coordinates
        px = (pts_nm[:, 0] - x_min) / nm_per_pixel
        py = (pts_nm[:, 1] - y_min) / nm_per_pixel
        _fill_polygon(mask, px, py)

    print(f"Loaded GDS: {filename} -> {width} x {height} pixels "
          f"({nm_per_pixel} nm/px, layer {layer})")
    return mask, (x_min, y_min)


def _fill_polygon(mask, px, py):
    """Scanline-fill a polygon into a 2D mask array."""
    from PIL import ImageDraw
    h, w = mask.shape
    img = Image.new('L', (w, h), 0)
    draw = ImageDraw.Draw(img)
    coords = list(zip(px.tolist(), py.tolist()))
    draw.polygon(coords, fill=255)
    arr = np.array(img, dtype=np.float64) / 255.0
    mask[:] = np.maximum(mask, arr)


def load_bitmap_mask(filename):
    """Load a bitmap mask (PBM, PNG, etc.) via Pillow.

    Returns array with 1.0 = feature (dark/chrome), 0.0 = clear.
    """
    img = Image.open(filename).convert('1')
    arr = np.array(img, dtype=np.uint8)
    # Pillow: black=0, white=1. PBM convention: 1=black=feature.
    mask = (1 - arr).astype(np.float64)
    h, w = mask.shape
    print(f"Loaded {w} x {h} bitmap mask: {filename}")
    return mask


def load_glp(filename, nm_per_pixel=1.0, padding=0):
    """Load an ICCAD 2013 .glp layout file.

    GLP files contain RECT and PGON entries with coordinates in nm.
    Format:
        RECT N M1  x y width height
        PGON N M1  x1 y1 x2 y2 ... xN yN
    """
    rects = []
    pgons = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('RECT'):
                parts = line.split()
                # RECT N M1 x y w h
                x, y, w, h = float(parts[3]), float(parts[4]), float(parts[5]), float(parts[6])
                rects.append((x, y, w, h))
            elif line.startswith('PGON'):
                parts = line.split()
                # PGON N M1 x1 y1 x2 y2 ...
                coords = [float(v) for v in parts[3:]]
                points = [(coords[i], coords[i+1]) for i in range(0, len(coords), 2)]
                pgons.append(points)

    # collect all points to find bounding box
    all_x, all_y = [], []
    for x, y, w, h in rects:
        all_x.extend([x, x + w])
        all_y.extend([y, y + h])
    for pts in pgons:
        for px, py in pts:
            all_x.append(px)
            all_y.append(py)

    x_min = min(all_x) - padding * nm_per_pixel
    y_min = min(all_y) - padding * nm_per_pixel
    x_max = max(all_x) + padding * nm_per_pixel
    y_max = max(all_y) + padding * nm_per_pixel

    width = int(np.ceil((x_max - x_min) / nm_per_pixel))
    height = int(np.ceil((y_max - y_min) / nm_per_pixel))
    mask = np.zeros((height, width), dtype=np.float64)

    # rasterize rectangles
    from PIL import ImageDraw
    img = Image.new('L', (width, height), 0)
    draw = ImageDraw.Draw(img)

    for x, y, w, h in rects:
        px1 = (x - x_min) / nm_per_pixel
        py1 = (y - y_min) / nm_per_pixel
        px2 = (x + w - x_min) / nm_per_pixel
        py2 = (y + h - y_min) / nm_per_pixel
        draw.rectangle([px1, py1, px2, py2], fill=255)

    for pts in pgons:
        coords = [((px - x_min) / nm_per_pixel, (py - y_min) / nm_per_pixel) for px, py in pts]
        draw.polygon(coords, fill=255)

    mask = np.array(img, dtype=np.float64) / 255.0
    print(f"Loaded GLP: {filename} -> {width} x {height} pixels "
          f"({nm_per_pixel} nm/px, {len(rects)} rects, {len(pgons)} pgons)")
    return mask


def load_mask(filename, layer=None, nm_per_pixel=1.0, padding=20):
    """Auto-detect format and load mask."""
    ext = os.path.splitext(filename)[1].lower()
    if ext in ('.gds', '.gdsii', '.gds2', '.oas'):
        mask, origin = load_gds(filename, layer=layer,
                                nm_per_pixel=nm_per_pixel, padding=padding)
        return mask
    elif ext == '.glp':
        return load_glp(filename, nm_per_pixel=nm_per_pixel, padding=padding)
    else:
        return load_bitmap_mask(filename)


# ======================================================================
# Illumination source shapes
# ======================================================================

def make_source_conventional(N, sigma):
    """Conventional (disk) partially-coherent source.

    Parameters
    ----------
    N : int
        Grid size (NxN in frequency space, should be odd).
    sigma : float
        Partial coherence factor (0 = coherent, 1 = edge of pupil).

    Returns
    -------
    source : ndarray (float64, 2D)
        Source intensity distribution in normalized frequency space.
        Coordinate range is [-1, 1] mapping to f/f_cutoff.
    """
    y, x = np.mgrid[0:N, 0:N]
    center = N // 2
    r = np.sqrt((x - center)**2 + (y - center)**2) / center
    return (r <= sigma).astype(np.float64)


def make_source_annular(N, sigma_inner, sigma_outer):
    """Annular illumination source."""
    y, x = np.mgrid[0:N, 0:N]
    center = N // 2
    r = np.sqrt((x - center)**2 + (y - center)**2) / center
    return ((r >= sigma_inner) & (r <= sigma_outer)).astype(np.float64)


def make_source_dipole(N, sigma, center_offset, radius, orientation='x'):
    """Dipole illumination (two poles along x or y axis).

    Parameters
    ----------
    N : int
        Grid size.
    sigma : float
        Not used directly; kept for API consistency.
    center_offset : float
        Distance of each pole center from origin (in normalized freq).
    radius : float
        Radius of each pole (in normalized freq).
    orientation : str
        'x' for horizontal dipole, 'y' for vertical.
    """
    y, x = np.mgrid[0:N, 0:N]
    c = N // 2
    fx = (x - c) / c
    fy = (y - c) / c

    if orientation == 'x':
        r1 = np.sqrt((fx - center_offset)**2 + fy**2)
        r2 = np.sqrt((fx + center_offset)**2 + fy**2)
    else:
        r1 = np.sqrt(fx**2 + (fy - center_offset)**2)
        r2 = np.sqrt(fx**2 + (fy + center_offset)**2)

    return ((r1 <= radius) | (r2 <= radius)).astype(np.float64)


def make_source_quadrupole(N, center_offset, radius):
    """Quadrupole illumination (four poles at 45-degree positions)."""
    y, x = np.mgrid[0:N, 0:N]
    c = N // 2
    fx = (x - c) / c
    fy = (y - c) / c

    d = center_offset / np.sqrt(2)
    source = np.zeros((N, N), dtype=np.float64)
    for sx, sy in [(d, d), (d, -d), (-d, d), (-d, -d)]:
        r = np.sqrt((fx - sx)**2 + (fy - sy)**2)
        source = np.maximum(source, (r <= radius).astype(np.float64))
    return source


# ======================================================================
# Pupil function and Zernike aberrations
# ======================================================================

def zernike(n, m, rho, theta):
    """Evaluate a single Zernike polynomial Z_n^m(rho, theta).

    Uses the standard (Noll) convention.  Only a handful of low-order
    terms are needed for typical lithography aberrations.
    """
    from scipy.special import factorial
    m_abs = abs(m)
    R = np.zeros_like(rho)
    for s in range(int((n - m_abs) / 2) + 1):
        num = (-1)**s * factorial(n - s)
        den = factorial(s) * factorial((n + m_abs) // 2 - s) * factorial((n - m_abs) // 2 - s)
        R += (num / den) * rho**(n - 2*s)

    if m >= 0:
        return R * np.cos(m_abs * theta)
    else:
        return R * np.sin(m_abs * theta)


def make_pupil(N, na, wavelength, aberrations=None):
    """Create the complex pupil function P(fx, fy).

    Parameters
    ----------
    N : int
        Grid size.
    na : float
        Numerical aperture.
    wavelength : float
        Wavelength in nm.
    aberrations : dict or None
        Zernike coefficients as {(n, m): coeff_in_waves, ...}.
        Example: {(2,0): 0.05} for 0.05 waves of defocus.

    Returns
    -------
    pupil : ndarray (complex128, 2D)
        Pupil function.  |P|=1 inside the aperture, 0 outside.
    """
    y, x = np.mgrid[0:N, 0:N]
    center = N // 2
    # normalized pupil coordinates: rho goes 0..1 at the edge of the aperture
    rho = np.sqrt((x - center)**2 + (y - center)**2) / center
    theta = np.arctan2(y - center, x - center)

    # circular aperture
    aperture = (rho <= 1.0).astype(np.float64)

    phase = np.zeros((N, N), dtype=np.float64)
    if aberrations:
        for (n, m), coeff in aberrations.items():
            phase += coeff * zernike(n, m, rho, theta) * 2 * np.pi  # convert waves to radians

    pupil = aperture * np.exp(1j * phase)
    return pupil


# ======================================================================
# Hopkins / SOCS optical model
# ======================================================================

def compute_tcc_kernels(source, pupil, n_kernels=None):
    """Compute SOCS kernels from the Hopkins TCC.

    The Transmission Cross-Coefficient is:
        TCC(f1, f2) = integral_source  P(s + f1) * P*(s + f2)  ds

    We decompose TCC via eigenvalue decomposition to get SOCS kernels
    (Sum of Coherent Systems):
        I(x) = sum_i  lambda_i * |phi_i * mask|^2

    For efficiency, we work in frequency domain on a discrete grid.

    Parameters
    ----------
    source : ndarray (float64, 2D)
        Source intensity distribution (NxN).
    pupil : ndarray (complex128, 2D)
        Pupil function (NxN).
    n_kernels : int or None
        Number of SOCS kernels to keep.  If None, keeps all with
        eigenvalue > 1e-3 * max_eigenvalue.

    Returns
    -------
    kernels : list of ndarray (complex128, 2D)
        SOCS kernels (spatial domain), sorted by eigenvalue magnitude.
    eigenvalues : ndarray (float64, 1D)
        Corresponding eigenvalues.
    """
    N = source.shape[0]
    assert source.shape == pupil.shape == (N, N)

    # Build TCC matrix.  For an NxN frequency grid, TCC is N^2 x N^2.
    # This can be large, so we keep N modest (e.g. 31-65).
    #
    # TCC[f1, f2] = sum_s  source(s) * pupil(s + f1) * conj(pupil(s + f2))
    #
    # We compute this via correlation: for each source point s,
    # the contribution is outer product of shifted-pupil vectors.

    # Flatten frequency indices
    N2 = N * N

    # Precompute pupil shifted by each source point using circular shifts.
    # For each source point (sy, sx) where source > 0, shift pupil.
    src_points = np.argwhere(source > 0)

    tcc = np.zeros((N2, N2), dtype=np.complex128)

    for sy, sx in src_points:
        # shift pupil by (sy - N//2, sx - N//2)
        dy = sy - N // 2
        dx = sx - N // 2
        shifted = np.roll(np.roll(pupil, -dy, axis=0), -dx, axis=1)
        svec = shifted.ravel()
        weight = source[sy, sx]
        tcc += weight * np.outer(svec, svec.conj())

    # normalize
    tcc /= source.sum()

    # Eigendecompose (TCC is Hermitian positive semi-definite)
    eigenvalues, eigenvectors = np.linalg.eigh(tcc)

    # Sort by descending eigenvalue
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx].real
    eigenvectors = eigenvectors[:, idx]

    # Truncate
    max_ev = eigenvalues[0]
    if n_kernels is None:
        keep = np.sum(eigenvalues > 1e-3 * max_ev)
        n_kernels = max(int(keep), 1)
    n_kernels = min(n_kernels, len(eigenvalues))

    kernels = []
    kept_evals = eigenvalues[:n_kernels]
    for i in range(n_kernels):
        k = eigenvectors[:, i].reshape(N, N)
        kernels.append(k)

    print(f"SOCS: keeping {n_kernels} kernels "
          f"(top eigenvalues: {kept_evals[:5]})")

    return kernels, kept_evals


def compute_aerial_image_socs(mask, kernels, eigenvalues):
    """Compute aerial image using SOCS decomposition.

    I(x,y) = sum_i  lambda_i * |kernel_i ** mask|^2

    where ** denotes 2D convolution.
    """
    aerial = np.zeros_like(mask, dtype=np.float64)
    for lam, kernel in zip(eigenvalues, kernels):
        if lam < 0:
            continue
        field = fftconvolve(mask, kernel, mode='same')
        aerial += lam * np.abs(field)**2
    return aerial


def compute_aerial_image_coherent(mask, kernel):
    """Simple fully-coherent model: I = |h * mask|^2 (original behavior)."""
    field = fftconvolve(mask, kernel, mode='same')
    return field * field


# ======================================================================
# Jinc kernel (for coherent / legacy mode)
# ======================================================================

def make_jinc_kernel(filter_size, scale):
    """Create a Jinc (Airy) kernel for fully-coherent Hopkins imaging."""
    P = filter_size if filter_size % 2 == 1 else filter_size + 1
    offset = P // 2
    y, x = np.mgrid[0:P, 0:P]
    r = np.sqrt((y - offset)**2 + (x - offset)**2).astype(np.float64)

    twopi = 2.0 * np.pi
    with np.errstate(divide='ignore', invalid='ignore'):
        kernel = bessj1(twopi * r * scale) / (twopi * r * scale)
    kernel[offset, offset] = 0.5
    kernel /= kernel.sum()
    return kernel


# ======================================================================
# Resist model
# ======================================================================

def apply_resist(aerial, threshold=0.2, blur_nm=0.0, nm_per_pixel=1.0):
    """Lumped-parameter resist model.

    Parameters
    ----------
    aerial : ndarray (float64, 2D)
        Aerial image intensity.
    threshold : float
        Resist clearing threshold (fraction of normalized intensity).
    blur_nm : float
        Acid diffusion length in nm (Gaussian sigma).  0 = no blur.
    nm_per_pixel : float
        Pixel size in nm (needed to convert blur_nm to pixels).

    Returns
    -------
    resist : ndarray (float64, 2D)
        Resist image after blur (continuous, 0..1).
    contour : ndarray (uint8, 2D)
        Binary contour (thresholded resist image).
    """
    if blur_nm > 0 and nm_per_pixel > 0:
        sigma_px = blur_nm / nm_per_pixel
        resist = gaussian_filter(aerial, sigma=sigma_px)
    else:
        resist = aerial.copy()

    contour = (resist > threshold).astype(np.uint8)
    return resist, contour


# ======================================================================
# Output helpers
# ======================================================================

def save_image(filename, data, normalize=True):
    """Save a 2D array as an image (format inferred from extension)."""
    if normalize:
        mn, mx = data.min(), data.max()
        if mx - mn > 0:
            arr = ((data - mn) / (mx - mn) * 255).astype(np.uint8)
        else:
            arr = np.zeros_like(data, dtype=np.uint8)
    else:
        arr = np.clip(data * 255, 0, 255).astype(np.uint8)
    Image.fromarray(arr, mode='L').save(filename)
    h, w = data.shape
    print(f"Saved {w} x {h}: {filename}")


def save_binary(filename, data):
    """Save binary mask (1 = feature = black in PBM)."""
    if data.dtype != np.uint8:
        data = (data > 0.5).astype(np.uint8)
    img = Image.fromarray((1 - data) * 255, mode='L').convert('1')
    img.save(filename)
    h, w = data.shape
    print(f"Saved {w} x {h}: {filename}")


# ======================================================================
# Mask to GDS (vectorization)
# ======================================================================

def _simplify_polygon(points, tolerance):
    """Simplify a polygon using the Ramer-Douglas-Peucker algorithm."""
    if len(points) < 3:
        return points

    # find the point farthest from the line between first and last
    start, end = points[0], points[-1]
    line_vec = end - start
    line_len = np.linalg.norm(line_vec)

    if line_len < 1e-10:
        dists = np.linalg.norm(points - start, axis=1)
    else:
        line_unit = line_vec / line_len
        proj = np.dot(points - start, line_unit)
        closest = start + np.outer(proj, line_unit)
        dists = np.linalg.norm(points - closest, axis=1)

    max_idx = np.argmax(dists)
    max_dist = dists[max_idx]

    if max_dist > tolerance:
        left = _simplify_polygon(points[:max_idx + 1], tolerance)
        right = _simplify_polygon(points[max_idx:], tolerance)
        return np.vstack([left[:-1], right])
    else:
        return np.array([start, end])


def _rectilinearize(points, snap_threshold=0.5):
    """Snap near-axis-aligned edges to be exactly rectilinear.

    For each edge, if the dx or dy component is small relative to the
    other, snap it to zero.  This produces Manhattan geometry preferred
    by mask shops.
    """
    result = points.copy()
    n = len(result)
    for i in range(n):
        j = (i + 1) % n
        dx = abs(result[j, 0] - result[i, 0])
        dy = abs(result[j, 1] - result[i, 1])
        edge_len = max(dx, dy, 1e-10)
        if dx / edge_len < snap_threshold:
            avg_x = (result[i, 0] + result[j, 0]) / 2
            result[i, 0] = avg_x
            result[j, 0] = avg_x
        if dy / edge_len < snap_threshold:
            avg_y = (result[i, 1] + result[j, 1]) / 2
            result[i, 1] = avg_y
            result[j, 1] = avg_y
    return result


def mask_to_polygons(mask, nm_per_pixel=1.0, origin=(0.0, 0.0),
                     simplify_tolerance=1.0, rectilinear=True):
    """Convert a binary mask to a list of polygon coordinate arrays.

    Uses marching squares to find contours, then simplifies them.

    Parameters
    ----------
    mask : ndarray (2D)
        Binary mask (>0.5 = feature).
    nm_per_pixel : float
        Pixel size in nm (for coordinate scaling).
    origin : tuple (x_nm, y_nm)
        World-coordinate offset of pixel (0,0).
    simplify_tolerance : float
        Douglas-Peucker simplification tolerance in pixels.
        Higher = fewer vertices = simpler shapes.
    rectilinear : bool
        If True, snap near-Manhattan edges to exact Manhattan.

    Returns
    -------
    polygons : list of ndarray, each shape (N, 2) in nm coordinates.
    """
    binary = (mask > 0.5).astype(np.float64)
    # pad to ensure closed contours at boundaries
    padded = np.pad(binary, 1, mode='constant', constant_values=0)
    contours = measure.find_contours(padded, 0.5)

    polygons = []
    for contour in contours:
        # remove padding offset, contour is (row, col) = (y, x)
        pts = contour - 1.0  # undo padding

        if simplify_tolerance > 0:
            pts = _simplify_polygon(pts, simplify_tolerance)

        if len(pts) < 3:
            continue

        if rectilinear:
            pts = _rectilinearize(pts)

        # convert from (row, col) pixel coords to (x_nm, y_nm)
        xy_nm = np.column_stack([
            pts[:, 1] * nm_per_pixel + origin[0],
            pts[:, 0] * nm_per_pixel + origin[1],
        ])
        polygons.append(xy_nm)

    return polygons


def save_gds(filename, polygons, layer=0, datatype=0, cell_name='TOP',
             units=(1e-9, 1e-12)):
    """Write polygons to a GDS file.

    Parameters
    ----------
    filename : str
        Output .gds path.
    polygons : list of ndarray
        Polygon vertices in nm coordinates.
    layer : int
        GDS layer number.
    datatype : int
        GDS datatype.
    cell_name : str
        Top cell name.
    units : tuple
        (user_unit, db_unit) in meters. Default: 1nm user, 1pm database.
    """
    import gdstk

    lib = gdstk.Library(unit=units[0], precision=units[1])
    cell = lib.new_cell(cell_name)

    for poly_nm in polygons:
        # gdstk uses library units; with unit=1e-9, coordinates are in nm
        cell.add(gdstk.Polygon(poly_nm, layer=layer, datatype=datatype))

    lib.write_gds(filename)
    print(f"Saved GDS: {filename} ({len(polygons)} polygons, layer {layer}/{datatype})")


def polygon_complexity(polygons):
    """Compute a shape complexity metric for a set of polygons.

    Returns a score where lower = simpler/more manufacturable:
        total_vertices + polygon_count * 4

    This penalizes both many polygons and complex polygon shapes.
    Rectangles (4 vertices) are cheapest.
    """
    if not polygons:
        return 0
    total_verts = sum(len(p) for p in polygons)
    n_polys = len(polygons)
    return total_verts + n_polys * 4


# ======================================================================
# Mask metrics (for OPC)
# ======================================================================

def mask_diff(a, b):
    return int(np.sum((a > 0.5) != (b > 0.5)))


def mask_complexity(mask):
    m = (mask > 0.5).astype(np.uint8)
    return int(np.sum(m[:, 1:] != m[:, :-1]) + np.sum(m[1:, :] != m[:-1, :]))


def neighbor_count_8(mask, x, y):
    h, w = mask.shape
    val = mask[y, x] > 0.5
    cnt = 0
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            ny = max(0, min(h - 1, y + dy))
            nx = max(0, min(w - 1, x + dx))
            if (mask[ny, nx] > 0.5) == val:
                cnt += 1
    return cnt


# ======================================================================
# OPC via simulated annealing
# ======================================================================

INITIAL_TEMP = 50000
FINAL_TEMP = 1.0
ERROR_WEIGHT = 1.0
COMPLEXITY_WEIGHT = 3.0
SHAPE_WEIGHT = 1.0  # polygon/vertex simplicity penalty
CLIMB_WEIGHT = 50.0


def _schedule(temp):
    if temp > 1000:
        return 0.95 * temp
    elif temp > 100:
        return 0.99 * temp
    else:
        return 0.98 * temp


def _compute_opc_cost(opced_mask, target, compute_aerial_fn, threshold,
                      initial_error, optimal_complexity, initial_shape_cost,
                      simplify_tolerance=2.0):
    """Compute the full OPC cost including pixel error, transition complexity,
    and polygon shape simplicity."""
    aerial = compute_aerial_fn(opced_mask)
    contour = (aerial > threshold).astype(np.float64)
    err = mask_diff(contour, target) / initial_error
    cplx = mask_complexity(opced_mask) / optimal_complexity

    # shape cost: evaluated less frequently (expensive), so cache-friendly
    # polygon complexity encourages fewer, simpler (rectangular) shapes
    polys = mask_to_polygons(opced_mask, simplify_tolerance=simplify_tolerance)
    shape = polygon_complexity(polys) / max(initial_shape_cost, 1)

    return (ERROR_WEIGHT * err +
            COMPLEXITY_WEIGHT * cplx +
            SHAPE_WEIGHT * shape)


def anneal(opced_mask, target, compute_aerial_fn, threshold=0.2,
           simplify_tolerance=2.0, rng=None):
    """Pixel-based OPC via simulated annealing.

    The cost function balances three terms:
    1. Aerial image error vs target (pixel accuracy)
    2. Pixel-level transition complexity (fewer transitions = smoother)
    3. Polygon shape simplicity (fewer polygons, fewer vertices)

    Parameters
    ----------
    opced_mask : ndarray (float64, 2D)
        Initial mask (modified in place).
    target : ndarray (float64, 2D)
        Desired output pattern.
    compute_aerial_fn : callable
        Function: mask -> aerial_image.
    threshold : float
        Resist threshold.
    simplify_tolerance : float
        Douglas-Peucker tolerance for polygon simplification (pixels).
    """
    if rng is None:
        rng = np.random.default_rng(2)

    h, w = target.shape

    aerial = compute_aerial_fn(opced_mask)
    contour = (aerial > threshold).astype(np.float64)

    initial_error = max(mask_diff(contour, target), 1)
    optimal_complexity = max(mask_complexity(opced_mask), 1)
    initial_polys = mask_to_polygons(opced_mask, simplify_tolerance=simplify_tolerance)
    initial_shape_cost = max(polygon_complexity(initial_polys), 1)

    old_cost = ERROR_WEIGHT + COMPLEXITY_WEIGHT + SHAPE_WEIGHT

    temperature = INITIAL_TEMP
    accept = climb = reject = 0
    counter = 20
    # full shape cost is expensive; compute it every N outer iterations
    shape_eval_interval = 5
    shape_counter = 0
    cached_shape_cost = 1.0  # normalized

    print(f"OPC: initial error={initial_error}  complexity={optimal_complexity}  "
          f"shape_cost={initial_shape_cost}")

    while temperature > FINAL_TEMP:
        coarseness = 4 if temperature > 5000 else 3 if temperature > 1000 else 2 if temperature > 100 else 1

        if counter == 20:
            cur_err = mask_diff((compute_aerial_fn(opced_mask) > threshold).astype(np.float64), target)
            cur_polys = mask_to_polygons(opced_mask, simplify_tolerance=simplify_tolerance)
            cur_shape = polygon_complexity(cur_polys)
            cached_shape_cost = cur_shape / max(initial_shape_cost, 1)
            print(f"\nT={temperature:.1e} sz={coarseness} "
                  f"acc={accept} clmb={climb} rej={reject} err={cur_err} "
                  f"shapes={len(cur_polys)} verts={sum(len(p) for p in cur_polys)} "
                  f"cost={old_cost:.3f}")
            counter = 0
            accept = climb = reject = 0
        else:
            counter += 1

        print(".", end="", flush=True)

        # periodically update shape cost
        shape_counter += 1
        if shape_counter >= shape_eval_interval:
            polys = mask_to_polygons(opced_mask, simplify_tolerance=simplify_tolerance)
            cached_shape_cost = polygon_complexity(polys) / max(initial_shape_cost, 1)
            shape_counter = 0

        for _ in range(1000):
            while True:
                x, y = rng.integers(0, w), rng.integers(0, h)
                nc = neighbor_count_8(opced_mask, x, y)
                if rng.random() >= (9.0 - nc * 0.1):
                    break

            x2, y2 = min(w, x + 3), min(h, y + 3)
            opced_mask[y:y2, x:x2] = 1.0 - opced_mask[y:y2, x:x2]

            new_aerial = compute_aerial_fn(opced_mask)
            new_contour = (new_aerial > threshold).astype(np.float64)
            new_cost = (ERROR_WEIGHT * mask_diff(new_contour, target) / initial_error +
                        COMPLEXITY_WEIGHT * mask_complexity(opced_mask) / optimal_complexity +
                        SHAPE_WEIGHT * cached_shape_cost)
            delta_c = new_cost - old_cost

            if delta_c < 0:
                accept += 1
                old_cost = new_cost
            elif rng.random() < np.exp(-CLIMB_WEIGHT * INITIAL_TEMP * delta_c / temperature):
                climb += 1
                old_cost = new_cost
            else:
                reject += 1
                opced_mask[y:y2, x:x2] = 1.0 - opced_mask[y:y2, x:x2]

        temperature = _schedule(temperature)

    final_contour = (compute_aerial_fn(opced_mask) > threshold).astype(np.float64)
    final_polys = mask_to_polygons(opced_mask, simplify_tolerance=simplify_tolerance)
    print(f"\nOPC done: error={mask_diff(final_contour, target)} "
          f"shapes={len(final_polys)} verts={sum(len(p) for p in final_polys)} "
          f"complexity={mask_complexity(opced_mask)}")
    return opced_mask


# ======================================================================
# High-level simulation API
# ======================================================================

class LithoSim:
    """Lithography simulation engine.

    Parameters
    ----------
    na : float
        Numerical aperture.
    wavelength : float
        Exposure wavelength in nm.
    nm_per_pixel : float
        Mask pixel size in nm.
    sigma : float
        Partial coherence factor (conventional source).
    source_type : str
        'conventional', 'annular', 'dipole_x', 'dipole_y', 'quadrupole', 'coherent'.
    sigma_inner : float
        Inner radius for annular source.
    sigma_outer : float
        Outer radius for annular source.
    n_kernels : int or None
        Number of SOCS kernels to retain.
    threshold : float
        Resist clearing threshold.
    resist_blur_nm : float
        Acid diffusion length (Gaussian sigma in nm).
    aberrations : dict or None
        Zernike aberration coefficients {(n,m): waves, ...}.
    defocus_nm : float
        Defocus in nm (converted to Zernike Z_2^0 coefficient).
    """

    def __init__(self, na=0.95, wavelength=193.0, nm_per_pixel=8.0,
                 sigma=0.75, source_type='conventional',
                 sigma_inner=0.5, sigma_outer=0.9,
                 n_kernels=None, threshold=0.2, resist_blur_nm=0.0,
                 aberrations=None, defocus_nm=0.0):
        self.na = na
        self.wavelength = wavelength
        self.nm_per_pixel = nm_per_pixel
        self.sigma = sigma
        self.source_type = source_type
        self.sigma_inner = sigma_inner
        self.sigma_outer = sigma_outer
        self.n_kernels = n_kernels
        self.threshold = threshold
        self.resist_blur_nm = resist_blur_nm
        self.aberrations = aberrations or {}
        self.defocus_nm = defocus_nm

        # derived
        self.min_feature_px = (wavelength / na) / nm_per_pixel
        self._kernels = None
        self._eigenvalues = None
        self._jinc_kernel = None

    def _setup_optics(self):
        """Build SOCS kernels or Jinc kernel depending on source type."""
        if self.source_type == 'coherent':
            scale = (self.na / self.wavelength) * self.nm_per_pixel
            self._jinc_kernel = make_jinc_kernel(int(self.min_feature_px), scale)
            print(f"Coherent mode: Jinc kernel {self._jinc_kernel.shape[0]}x{self._jinc_kernel.shape[1]}")
            return

        # frequency grid size — must be large enough to capture the kernel
        # but small enough for tractable TCC eigendecomposition
        N = int(self.min_feature_px)
        N = max(N, 15)
        if N % 2 == 0:
            N += 1
        # cap at reasonable size for eigendecomposition (N^2 x N^2 matrix)
        N = min(N, 51)

        print(f"Building SOCS model: N={N}, source={self.source_type}, sigma={self.sigma}")

        # build illumination source
        if self.source_type == 'conventional':
            source = make_source_conventional(N, self.sigma)
        elif self.source_type == 'annular':
            source = make_source_annular(N, self.sigma_inner, self.sigma_outer)
        elif self.source_type == 'dipole_x':
            source = make_source_dipole(N, self.sigma, 0.7, 0.2, 'x')
        elif self.source_type == 'dipole_y':
            source = make_source_dipole(N, self.sigma, 0.7, 0.2, 'y')
        elif self.source_type == 'quadrupole':
            source = make_source_quadrupole(N, 0.7, 0.2)
        else:
            raise ValueError(f"Unknown source type: {self.source_type}")

        # build pupil with aberrations
        ab = dict(self.aberrations)
        if self.defocus_nm != 0:
            # defocus -> Zernike Z(2,0) coefficient in waves
            # W_defocus = defocus * NA^2 / (2 * wavelength)  (paraxial approx)
            defocus_waves = self.defocus_nm * self.na**2 / (2 * self.wavelength)
            ab[(2, 0)] = ab.get((2, 0), 0) + defocus_waves

        pupil = make_pupil(N, self.na, self.wavelength, ab if ab else None)

        self._kernels, self._eigenvalues = compute_tcc_kernels(
            source, pupil, n_kernels=self.n_kernels)

    def aerial_image(self, mask):
        """Compute aerial image for the given binary mask."""
        if self._kernels is None and self._jinc_kernel is None:
            self._setup_optics()

        if self.source_type == 'coherent':
            return compute_aerial_image_coherent(mask, self._jinc_kernel)
        else:
            return compute_aerial_image_socs(mask, self._kernels, self._eigenvalues)

    def resist_image(self, aerial):
        """Apply resist model to aerial image."""
        return apply_resist(aerial, self.threshold, self.resist_blur_nm,
                            self.nm_per_pixel)

    def simulate(self, mask):
        """Full simulation: mask -> aerial -> resist -> contour.

        Returns (aerial, resist, contour).
        """
        aerial = self.aerial_image(mask)
        resist, contour = self.resist_image(aerial)
        return aerial, resist, contour

    def opc(self, target_mask):
        """Run OPC on the target mask. Returns optimized mask."""
        if self._kernels is None and self._jinc_kernel is None:
            self._setup_optics()
        opced = target_mask.copy()
        return anneal(opced, target_mask, self.aerial_image, self.threshold)


# ======================================================================
# CLI
# ======================================================================

def parse_layer(s):
    """Parse layer spec like '1/0' or '5' into (layer, datatype) tuple."""
    parts = s.split('/')
    layer = int(parts[0])
    datatype = int(parts[1]) if len(parts) > 1 else 0
    return (layer, datatype)


def main():
    parser = argparse.ArgumentParser(
        description='Lithography simulator with GDS support and Hopkins/SOCS optics.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s mask.pbm out                          # coherent sim, bitmap input
  %(prog)s --sigma 0.7 layout.gds out            # partial coherence, GDS input
  %(prog)s --layer 1/0 --nm-per-pixel 4 chip.gds out
  %(prog)s --source annular --sigma-inner 0.5 --sigma-outer 0.9 chip.gds out
  %(prog)s --defocus 50 --resist-blur 30 chip.gds out
  %(prog)s --opc mask.pbm out                    # run OPC
""")

    parser.add_argument('input', help='Input mask file (GDS, OASIS, PBM, PNG)')
    parser.add_argument('output_prefix', help='Output prefix (e.g. results/sim)')

    # Optics
    parser.add_argument('--na', type=float, default=0.95, help='Numerical aperture (default: 0.95)')
    parser.add_argument('--wavelength', type=float, default=193.0, help='Wavelength in nm (default: 193)')
    parser.add_argument('--nm-per-pixel', type=float, default=8.0, help='Resolution in nm/pixel (default: 8)')
    parser.add_argument('--source', default='conventional',
                        choices=['coherent', 'conventional', 'annular', 'dipole_x', 'dipole_y', 'quadrupole'],
                        help='Illumination source type (default: conventional)')
    parser.add_argument('--sigma', type=float, default=0.75, help='Partial coherence sigma (default: 0.75)')
    parser.add_argument('--sigma-inner', type=float, default=0.5, help='Annular inner sigma')
    parser.add_argument('--sigma-outer', type=float, default=0.9, help='Annular outer sigma')
    parser.add_argument('--n-kernels', type=int, default=None, help='Number of SOCS kernels (auto if not set)')
    parser.add_argument('--defocus', type=float, default=0.0, help='Defocus in nm')

    # Resist
    parser.add_argument('--threshold', type=float, default=0.2, help='Resist threshold (default: 0.2)')
    parser.add_argument('--resist-blur', type=float, default=0.0, help='Resist acid diffusion blur in nm')

    # GDS options
    parser.add_argument('--layer', type=str, default=None, help='GDS layer/datatype (e.g. 1/0)')
    parser.add_argument('--padding', type=int, default=20, help='Padding around GDS bbox in pixels')

    # OPC
    parser.add_argument('--opc', action='store_true', help='Run OPC optimization')

    args = parser.parse_args()

    t0 = time.time()

    layer = parse_layer(args.layer) if args.layer else None
    mask = load_mask(args.input, layer=layer,
                     nm_per_pixel=args.nm_per_pixel, padding=args.padding)

    sim = LithoSim(
        na=args.na,
        wavelength=args.wavelength,
        nm_per_pixel=args.nm_per_pixel,
        sigma=args.sigma,
        source_type=args.source,
        sigma_inner=args.sigma_inner,
        sigma_outer=args.sigma_outer,
        n_kernels=args.n_kernels,
        threshold=args.threshold,
        resist_blur_nm=args.resist_blur,
        defocus_nm=args.defocus,
    )

    os.makedirs(os.path.dirname(args.output_prefix) or '.', exist_ok=True)

    if args.opc:
        opced_mask = sim.opc(mask)
        save_binary(f"{args.output_prefix}_opc_mask.pbm", (opced_mask > 0.5).astype(np.uint8))
        aerial, resist, contour = sim.simulate(opced_mask)
        # export OPC mask as GDS
        polys = mask_to_polygons(opced_mask, nm_per_pixel=args.nm_per_pixel,
                                 simplify_tolerance=2.0, rectilinear=True)
        out_layer = parse_layer(args.layer) if args.layer else (0, 0)
        save_gds(f"{args.output_prefix}_opc_mask.gds", polys,
                 layer=out_layer[0], datatype=out_layer[1])
    else:
        aerial, resist, contour = sim.simulate(mask)

    save_image(f"{args.output_prefix}_aerial.png", aerial)
    save_image(f"{args.output_prefix}_resist.png", resist)
    save_binary(f"{args.output_prefix}_contour.pbm", contour)

    # also export contour as GDS
    contour_polys = mask_to_polygons(contour.astype(np.float64),
                                     nm_per_pixel=args.nm_per_pixel,
                                     simplify_tolerance=1.0, rectilinear=True)
    out_layer = parse_layer(args.layer) if args.layer else (0, 0)
    save_gds(f"{args.output_prefix}_contour.gds", contour_polys,
             layer=out_layer[0], datatype=out_layer[1])

    print(f"Total: {time.time()-t0:.2f}s")


if __name__ == "__main__":
    main()
