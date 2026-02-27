from streamline_extraction import StrealimeExtractor
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.optimize import minimize

import torch


def get_streamlines(extracted_data):
    streamlines = {}

    for surface_tag in extracted_data['surfaces'].keys():
        surface_data = extracted_data['surfaces'][surface_tag]

        streamlines_surfaces_points = []
        streamlines_surfaces_splines = []

        for separatrix_id in surface_data['separatrices_ids']:

            points = extracted_data['vertices'][separatrix_id]
            # tck         = points_to_spline(points)

            streamlines_surfaces_points.append(points)
            # streamlines_surfaces_splines.append(tck)

        streamlines[surface_tag] = {
            'points': streamlines_surfaces_points,
            # 'splines' : streamlines_surfaces_splines,
        }

    return streamlines


streamline_extractor = StrealimeExtractor()
extracted_data = streamline_extractor.extract_mesh_from_msh(
    path_msh_file="./msh_files/remeshing.msh")
# extracted_data = streamline_extractor.extract_mesh_from_msh(path_msh_file="./msh_files/remeshing_v1.msh")

streamlines = get_streamlines(extracted_data)
streamlines[2].keys()

for key in streamlines.keys():
    for streamline in streamlines[key]['points']:
        # print(len(streamline))
        print((streamline.shape))


intersections = find_all_intersections()


def points_to_spline(points: np.ndarray, smoothing: float = 0):
    points = np.asarray(points, dtype=float)
    mask = np.concatenate(
        [[True], np.any(np.diff(points, axis=0) != 0, axis=1)])
    points = points[mask]
    if len(points) < 4:
        raise ValueError(f"Need at least 4 unique points, got {len(points)}")
    k = min(3, len(points) - 1)
    tck, _ = splprep([points[:, 0], points[:, 1],
                     points[:, 2]], s=smoothing, k=k)
    return tck


def _bboxes_gpu(tck, t_vals, device):
    pts = np.stack(splev(t_vals, tck), axis=1)
    pts = torch.tensor(pts, dtype=torch.float32, device=device)
    return torch.minimum(pts[:-1], pts[1:]), torch.maximum(pts[:-1], pts[1:])


def _bbox_candidates(all_tck: list, n: int, device: torch.device):
    t_vals = np.linspace(0, 1, n + 1)
    bboxes = [_bboxes_gpu(tck, t_vals, device) for tck in all_tck]
    candidates = {}

    for i in range(len(all_tck)):
        for j in range(i + 1, len(all_tck)):
            b1_min, b1_max = bboxes[i]
            b2_min, b2_max = bboxes[j]

            overlap = torch.all(
                (b1_min.unsqueeze(1) <= b2_max.unsqueeze(0)) &
                (b1_max.unsqueeze(1) >= b2_min.unsqueeze(0)),
                dim=2
            )

            ii, jj = torch.where(overlap)
            ii, jj = ii.cpu().numpy(), jj.cpu().numpy()
            if len(ii) > 0:
                candidates[(i, j)] = [
                    (t_vals[a], t_vals[a+1], t_vals[b], t_vals[b+1])
                    for a, b in zip(ii, jj)
                ]

    return candidates


def _refine(tck1, tck2, t1s, t1e, t2s, t2e, tol):
    def distance(params):
        p1 = np.array(splev(params[0], tck1))
        p2 = np.array(splev(params[1], tck2))
        return np.sum((p1 - p2) ** 2)

    result = minimize(
        distance,
        x0=[(t1s + t1e) / 2, (t2s + t2e) / 2],
        method='L-BFGS-B',
        bounds=[(t1s, t1e), (t2s, t2e)],
    )

    if result.fun < tol:
        t1_opt, t2_opt = result.x
        return {
            'point': np.array(splev(t1_opt, tck1)),
            't1': t1_opt,
            't2': t2_opt,
            'distance': np.sqrt(result.fun),
        }
    return None


def find_all_intersections(curves_points: list, n_bbox: int = 100,
                           tolerance: float = 1e-6, smoothing: float = 0):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    all_tck = [points_to_spline(pts, smoothing) for pts in curves_points]
    candidates = _bbox_candidates(all_tck, n_bbox, device)

    intersections = []
    for (i, j), segments in candidates.items():
        seen = []
        for (t1s, t1e, t2s, t2e) in segments:
            hit = _refine(all_tck[i], all_tck[j], t1s,
                          t1e, t2s, t2e, tolerance)
            if hit is None:
                continue
            if not any(np.linalg.norm(hit['point'] - p) < tolerance * 10 for p in seen):
                seen.append(hit['point'])
                intersections.append({'curve_i': i, 'curve_j': j, **hit})

    return intersections
