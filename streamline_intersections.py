import numpy as np
from scipy.interpolate import splprep, splev
from scipy.optimize import minimize

import torch


class StreamlineIntersections:

    def get_intersections(self, streamlines_all_surfaces):
        results = {}
        for surface_tag, data in streamlines_all_surfaces.items():
            streamlines = data['points']
            intersections = self.find_all_intersections(streamlines)
            split = self.split_streamlines_at_intersections(
                streamlines, intersections)
            results[surface_tag] = {
                'intersections': intersections,
                'streamlines_split': split
            }
        return results

    def points_to_spline(self, points: np.ndarray, smoothing: float = 0):

        points = np.asarray(points, dtype=float)
        mask = np.concatenate(
            [[True], np.any(np.diff(points, axis=0) != 0, axis=1)])
        points = points[mask]
        k = min(3, len(points) - 1)
        tck, _ = splprep([points[:, 0], points[:, 1],
                          points[:, 2]], s=smoothing, k=k)
        return tck

    def _bboxes_gpu(self, tck, t_vals, device):
        pts = np.stack(splev(t_vals, tck), axis=1)
        pts = torch.tensor(pts, dtype=torch.float32, device=device)
        return torch.minimum(pts[:-1], pts[1:]), torch.maximum(pts[:-1], pts[1:])

    def _bbox_candidates(self, all_tck: list, n: int, device: torch.device):
        t_vals = np.linspace(0, 1, n + 1)
        bboxes = [self._bboxes_gpu(tck, t_vals, device) for tck in all_tck]
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

    def split_streamlines_at_intersections(self, streamlines: list, intersections: list,
                                           n_samples_per_unit: int = 200) -> list:
        from collections import defaultdict
        from scipy.interpolate import splev

        split_ts = defaultdict(list)
        for hit in intersections:
            split_ts[hit['curve_i']].append(hit['t1'])
            split_ts[hit['curve_j']].append(hit['t2'])

        new_streamlines = []

        for idx, pts in enumerate(streamlines):
            tck = self.points_to_spline(pts, smoothing=0)
            k = tck[2]          # spline degree
            knots = tck[0]
            t_min = knots[k]            # valid range start
            t_max = knots[-(k + 1)]     # valid range end

            if idx not in split_ts:
                new_streamlines.append(pts)
                continue

            # Sort + deduplicate, clamp to valid range
            ts = sorted(set(
                np.clip(t, t_min, t_max) for t in split_ts[idx]
            ))

            boundaries = [t_min] + ts + [t_max]

            for k_seg in range(len(boundaries) - 1):
                t_start = boundaries[k_seg]
                t_end = boundaries[k_seg + 1]

                if t_end - t_start < 1e-9:
                    continue

                n = max(2, int((t_end - t_start) * n_samples_per_unit))
                t_samples = np.linspace(t_start, t_end, n)
                coords = np.array(splev(t_samples, tck)).T  # shape (n, D)

                new_streamlines.append(coords)

        return new_streamlines

    def _refine(self, tck1, tck2, t1s, t1e, t2s, t2e, tol):
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

    def find_all_intersections(self, curves_points: list, n_bbox: int = 100,
                               tolerance: float = 1e-6, smoothing: float = 0):
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        all_tck = [self.points_to_spline(pts, smoothing)
                   for pts in curves_points]
        candidates = self._bbox_candidates(all_tck, n_bbox, device)

        intersections = []
        for (i, j), segments in candidates.items():
            seen = []
            for (t1s, t1e, t2s, t2e) in segments:
                hit = self._refine(all_tck[i], all_tck[j], t1s,
                                   t1e, t2s, t2e, tolerance)
                if hit is None:
                    continue
                if not any(np.linalg.norm(hit['point'] - p) < tolerance * 10 for p in seen):
                    seen.append(hit['point'])
                    intersections.append({'curve_i': i, 'curve_j': j, **hit})

        return intersections
