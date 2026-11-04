"""Multi-root scan and continuation.

Requires the supplied stability.py and Function_determinant_abs_nov_EE.py.
The active MATLAB configuration is preserved, including bounds [-50, 80]
and [1, 50] Hz (the original introductory comment had outdated bounds).
"""
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.io import savemat
from scipy.optimize import minimize
from stability import Set_rates_IC
from Set_parameters_spatial import Set_parameters
from Function_determinant_abs_nov_EE import Function_determinant_abs_nov_EE


@dataclass
class Root:
    x: np.ndarray
    x_prev: object
    cost: float


def default_config():
    return SimpleNamespace(
        det_tol=1e-2, threshold_minimum=np.log(1e-2),
        global_re_bounds=np.array([-50., 80.]),
        global_im_bounds=2*np.pi*np.array([1., 50.]),
        max_active_roots=3,
        max_seed_attempts_initial=8, max_seed_attempts_audit=4,
        seed_re_half_width=5., seed_im_half_width=2*np.pi*5,
        cont_re_half_width=2.5, cont_im_half_width=2*np.pi*2.5,
        retry_re_half_width=5., retry_im_half_width=2*np.pi*5,
        merge_re_tol=0.5, merge_im_tol=2*np.pi*0.5,
        options_seed=dict(MaxIter=300, MaxFunEvals=500, TolX=1e-5,
                          TolFun=1e-6, Display='off'),
        options_continue=dict(MaxIter=120, MaxFunEvals=220, TolX=1e-5,
                              TolFun=1e-6, Display='off'),
    )


def empty_root_struct():
    return []


def clamp_point(x, re_bounds, im_bounds):
    return np.clip(np.asarray(x, float),
                   [re_bounds[0], im_bounds[0]],
                   [re_bounds[1], im_bounds[1]])


def bounded_cost(x, cost_func, re_bounds, im_bounds):
    if not (re_bounds[0] <= x[0] <= re_bounds[1]
            and im_bounds[0] <= x[1] <= im_bounds[1]):
        re_scale = max(re_bounds[1]-re_bounds[0], np.finfo(float).eps)
        im_scale = max(im_bounds[1]-im_bounds[0], np.finfo(float).eps)
        re_distance = max(re_bounds[0]-x[0], 0) + max(x[0]-re_bounds[1], 0)
        im_distance = max(im_bounds[0]-x[1], 0) + max(x[1]-im_bounds[1], 0)
        return 1e6 + (re_distance/re_scale)**2 + (im_distance/im_scale)**2
    return cost_func(x)


def refine_in_box(cost_func, starting_point, re_bounds, im_bounds, options):
    start = clamp_point(starting_point, re_bounds, im_bounds)
    result = minimize(
        bounded_cost, start, args=(cost_func, re_bounds, im_bounds),
        method='Nelder-Mead',
        options=dict(maxiter=options['MaxIter'], maxfev=options['MaxFunEvals'],
                     xatol=options['TolX'], fatol=options['TolFun'],
                     disp=False, adaptive=False))
    # SciPy status=0 means success; MATLAB's exitflag>0 means success.
    exitflag = 1 if result.success else 0
    return result.x, float(result.fun), exitflag


def is_valid_solution(solution, fval, exitflag, re_bounds, im_bounds,
                      threshold_minimum):
    return bool(exitflag > 0 and np.all(np.isfinite(solution))
                and np.isfinite(fval) and fval < threshold_minimum
                and re_bounds[0] <= solution[0] <= re_bounds[1]
                and im_bounds[0] <= solution[1] <= im_bounds[1])


def is_separated_seed(candidate, existing_seeds, separation_re, separation_im):
    if not len(existing_seeds):
        return True
    existing = np.asarray(existing_seeds)
    close = ((np.abs(existing[:, 0]-candidate[0]) < separation_re)
             & (np.abs(existing[:, 1]-candidate[1]) < separation_im))
    return not np.any(close)


def select_grid_seeds(grid_values, re_grid, im_grid, maximum_number,
                      separation_re, separation_im):
    n_im, n_re = grid_values.shape
    mask = np.zeros_like(grid_values, dtype=bool)
    for iim in range(n_im):
        for ire in range(n_re):
            neighborhood = grid_values[max(0, iim-1):min(n_im, iim+2),
                                       max(0, ire-1):min(n_re, ire+2)]
            mask[iim, ire] = grid_values[iim, ire] <= np.min(neighborhood)
    # MATLAB find/(:) use column-major order. Stable sorts preserve ties.
    flat = grid_values.ravel(order='F')
    local = np.flatnonzero(mask.ravel(order='F'))
    local = local[np.argsort(flat[local], kind='stable')]
    global_order = np.argsort(flat, kind='stable')
    seeds = []
    for candidates in (local, global_order):
        for index in candidates:
            iim, ire = np.unravel_index(index, grid_values.shape, order='F')
            candidate = [re_grid[ire], im_grid[iim], grid_values[iim, ire]]
            if is_separated_seed(candidate, seeds, separation_re, separation_im):
                seeds.append(candidate)
            if len(seeds) >= maximum_number:
                return np.asarray(seeds, float).reshape(-1, 3)
    return np.asarray(seeds, float).reshape(-1, 3)


def merge_root_sets(first_set, second_set, cfg):
    combined = deepcopy(first_set)
    if isinstance(second_set, Root):
        second_set = [second_set]
    for candidate in second_set:
        duplicate = next((j for j, root in enumerate(combined)
                          if abs(root.x[0]-candidate.x[0]) < cfg.merge_re_tol
                          and abs(root.x[1]-candidate.x[1]) < cfg.merge_im_tol), None)
        if duplicate is None:
            combined.append(deepcopy(candidate))
        else:
            existing_previous = combined[duplicate].x_prev
            if candidate.cost < combined[duplicate].cost:
                combined[duplicate] = deepcopy(candidate)
                if existing_previous is not None:
                    combined[duplicate].x_prev = existing_previous.copy()
            elif existing_previous is None and candidate.x_prev is not None:
                combined[duplicate].x_prev = candidate.x_prev.copy()
    return combined


def sort_and_trim_roots(roots, maximum_number):
    return sorted(roots, key=lambda root: root.x[0], reverse=True)[:maximum_number]


def roots_to_matrix(roots):
    return np.asarray([[root.x[0], root.x[1]/(2*np.pi), root.cost]
                       for root in roots], float).reshape(-1, 3)


def _window(center, re_half_width, im_half_width, cfg):
    return (np.array([max(cfg.global_re_bounds[0], center[0]-re_half_width),
                      min(cfg.global_re_bounds[1], center[0]+re_half_width)]),
            np.array([max(cfg.global_im_bounds[0], center[1]-im_half_width),
                      min(cfg.global_im_bounds[1], center[1]+im_half_width)]))


def scan_and_refine_roots(cost_func, re_grid, im_grid, max_seed_attempts, cfg):
    grid_values = np.empty((len(im_grid), len(re_grid)))
    for iim, im in enumerate(im_grid):
        for ire, re in enumerate(re_grid):
            grid_values[iim, ire] = cost_func(np.array([re, im]))
    best_cost_observed = float(np.min(grid_values))
    sep_re = 2*abs(re_grid[1]-re_grid[0]) if len(re_grid) > 1 else cfg.merge_re_tol
    sep_im = 2*abs(im_grid[1]-im_grid[0]) if len(im_grid) > 1 else cfg.merge_im_tol
    seeds = select_grid_seeds(grid_values, re_grid, im_grid,
                              max_seed_attempts, sep_re, sep_im)
    roots = empty_root_struct()
    for seed in seeds[:, :2]:
        rb, ib = _window(seed, cfg.seed_re_half_width, cfg.seed_im_half_width, cfg)
        solution, fval, flag = refine_in_box(cost_func, seed, rb, ib, cfg.options_seed)
        best_cost_observed = min(best_cost_observed, fval)
        if is_valid_solution(solution, fval, flag, rb, ib, cfg.threshold_minimum):
            roots = merge_root_sets(roots, Root(solution, None, fval), cfg)
    return sort_and_trim_roots(roots, cfg.max_active_roots), best_cost_observed


def continue_roots(cost_func, active_roots, cfg):
    continued = empty_root_struct()
    best_cost_observed = np.inf
    for root in active_roots:
        previous = root.x
        predictor = previous if root.x_prev is None else previous+(previous-root.x_prev)
        rb, ib = _window(previous, cfg.cont_re_half_width, cfg.cont_im_half_width, cfg)
        solution, fval, flag = refine_in_box(cost_func, predictor, rb, ib, cfg.options_continue)
        best_cost_observed = min(best_cost_observed, fval)
        valid = is_valid_solution(solution, fval, flag, rb, ib, cfg.threshold_minimum)
        if not valid:
            rb, ib = _window(previous, cfg.retry_re_half_width, cfg.retry_im_half_width, cfg)
            solution, fval, flag = refine_in_box(cost_func, previous, rb, ib, cfg.options_continue)
            best_cost_observed = min(best_cost_observed, fval)
            valid = is_valid_solution(solution, fval, flag, rb, ib, cfg.threshold_minimum)
        if valid:
            continued = merge_root_sets(continued, Root(solution, previous.copy(), fval), cfg)
    return sort_and_trim_roots(continued, cfg.max_active_roots), best_cost_observed


def _first_scalar(value, name):
    values = np.asarray(value, float).ravel()
    if not values.size:
        raise ValueError(f'{name} must not be empty')
    return float(values[0])


def Wrap_peak_anisotropy_EE(
        rho_EE_input, nu_from_input, g_from_input, v_fom_input, *,
        rho_II_input=0.0,
        rho_IE_input=0.0,
        rho_EI_input=0.0,
        param=None, output_dir='Results_5d_test', progress=True,
        n_u=21, n_v=21, n_initial_grid=50, n_audit_grid=18,
        u_max=1.5, v_max=1.5, audit_v_indices=None):
        
    """Run the original 21 x 21 spatial grid and save one .mat file per u.

    nu_from_input is the absolute external rate, not a threshold-rate ratio.
    v_fom_input (source spelling) is propagation velocity.

    Defaults are in
    Set_parameters_spatial.py. Optionally pass a complete param namespace.
    Lengths are mm, velocity is mm/s, and time is seconds.
    Caller parameters are copied before modification.

    n_u and n_v set the number of wavevector points from 0 to u_max/v_max.
    n_initial_grid is the resolution of the initial determinant scan at each u;
    n_audit_grid is the resolution of subsequent global audits. Both scans are
    square. audit_v_indices uses zero-based Python indices; by default it audits
    the midpoint and last v values, matching MATLAB indices 11 and 21 when
    n_v=21.

    Returns a list of the per-u dictionaries, including output_path.
    No computation is started merely by importing this module.
    """



    if n_u < 2 or n_v < 1 or n_initial_grid < 2 or n_audit_grid < 2:
        raise ValueError(
            "n_u, n_initial_grid and n_audit_grid must each be at least 2. "
            "n_v must be at least 1."
        )

    param = Set_parameters() if param is None else deepcopy(param)

    
    param.rho_EE = float(rho_EE_input)
    param.rho_II = float(rho_II_input)  
    param.rho_IE = float(rho_IE_input)
    param.rho_EI = float(rho_EI_input)


    rhos = (param.rho_EE, param.rho_II, param.rho_IE, param.rho_EI)
    if not all(-1.0 < rho < 1.0 for rho in rhos):
        raise ValueError("Each anisotropy parameter must lie strictly between -1 and 1.")
    
    param.Lbox = 20*param.l_ee
    List_u = np.linspace(0, float(u_max), int(n_u))
    List_v = np.linspace(0, float(v_max), int(n_v))
    param.velocity = _first_scalar(v_fom_input, 'v_fom_input')
    param.nu_ext = float(nu_from_input)
    param.gi = _first_scalar(g_from_input, 'g_from_input')
    param.ge = 1.07*param.gi
    # Store velocity explicitly as a float in the directory name; MATLAB's
    # %.d format in the original is ambiguous for noninteger velocities.


    folder = Path(output_dir)/(
    f'rhoEE{param.rho_EE:.1f}_rhoII{param.rho_II:.1f}_rhoIE{param.rho_IE:.1f}_rhoEI{param.rho_EI:.1f}_g{param.gi:.2f}_'
    f'nuext{param.nu_ext:.2f}_v{param.velocity:g}')

    folder.mkdir(parents=True, exist_ok=True)
    initial_nu = np.array([0.3, 1.])
    cfg = default_config()
    re_initial = np.linspace(*cfg.global_re_bounds, int(n_initial_grid))
    im_initial = np.linspace(*cfg.global_im_bounds, int(n_initial_grid))
    re_audit = np.linspace(*cfg.global_re_bounds, int(n_audit_grid))
    im_audit = np.linspace(*cfg.global_im_bounds, int(n_audit_grid))
    audit_indices = ({(int(n_v)-1)//2, int(n_v)-1}
                     if audit_v_indices is None
                     else {int(index) for index in audit_v_indices})
    if any(index < 0 or index >= n_v for index in audit_indices):
        raise ValueError("Each audit_v_indices value must be a valid zero-based v index")
    results = []
    for iu, u_param in enumerate(List_u):
        rates = Set_rates_IC(param, initial_nu)
        initial_nu = np.array([rates.nu_E, rates.nu_I])
        local_minima_all = np.empty((1, len(List_v)), dtype=object)
        all_roots = np.empty((1, len(List_v)), dtype=object)
        dominant_success = np.zeros((1, len(List_v)), dtype=bool)
        number_roots_found = np.zeros((1, len(List_v)), dtype=float)
        audit_was_run = np.zeros((1, len(List_v)), dtype=bool)
        active_roots = empty_root_struct()
        for iv, v_param in enumerate(List_v):
            param.u, param.v = u_param, v_param

            def cost_func(x):
                return Function_determinant_abs_nov_EE(complex(x[0], x[1]), param, rates)

            best_cost_observed = np.inf
            if iv == 0:
                active_roots, scan_best = scan_and_refine_roots(
                    cost_func, re_initial, im_initial, cfg.max_seed_attempts_initial, cfg)
                best_cost_observed = min(best_cost_observed, scan_best)
            else:
                continued, continue_best = continue_roots(cost_func, active_roots, cfg)
                best_cost_observed = min(best_cost_observed, continue_best)
                if iv in audit_indices or not continued:
                    audit_was_run[0, iv] = True
                    audited, audit_best = scan_and_refine_roots(
                        cost_func, re_audit, im_audit, cfg.max_seed_attempts_audit, cfg)
                    best_cost_observed = min(best_cost_observed, audit_best)
                    active_roots = merge_root_sets(continued, audited, cfg)
                else:
                    active_roots = continued
                active_roots = sort_and_trim_roots(active_roots, cfg.max_active_roots)
            active_roots = sort_and_trim_roots(active_roots, cfg.max_active_roots)
            matrix = roots_to_matrix(active_roots)
            all_roots[0, iv] = matrix
            number_roots_found[0, iv] = len(active_roots)
            if active_roots:
                local_minima_all[0, iv] = matrix[:1].copy()
                dominant_success[0, iv] = True
            else:
                if np.isinf(best_cost_observed):
                    best_cost_observed = np.nan
                local_minima_all[0, iv] = np.array([[np.nan, np.nan, best_cost_observed]])
            if progress:
                print(f'u {iu+1}/{len(List_u)}, v {iv+1}/{len(List_v)}: '
                      f'{len(active_roots)} retained roots', flush=True)
        path = folder/(f'Stability_rhoEE{param.rho_EE:.1f}_rhoII{param.rho_II:.1f}_rhoIE{param.rho_IE:.1f}_rhoEI{param.rho_EI:.1f}_nu_ext{param.nu_ext:.2f}_'
                       f'g{param.gi:.2f}_v{param.velocity:.1f}_kx{u_param:.2f}.mat')
        data = dict(List_v=List_v.reshape(1, -1), List_u=List_u.reshape(1, -1),
                    local_minima_all=local_minima_all, all_roots=all_roots,
                    dominant_success=dominant_success, number_roots_found=number_roots_found,
                    audit_was_run=audit_was_run, search_config=vars(cfg))
        savemat(path, data, long_field_names=True)
        results.append(dict(data, output_path=path))
    return results
