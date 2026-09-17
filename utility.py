import numpy as np
from scipy.signal import find_peaks


def oscillation_onset(
    time,
    values,
    prominence,
    min_cycles=3,
    period_tolerance=0.25,
    *,
    min_extrema=None,
):
    """Estimate the onset time of sustained oscillations.

    Sustained oscillations are identified from a sequence of
    alternating maxima and minima for which the periods between
    successive maxima and successive minima are approximately
    constant.

    The onset time is defined as one half of the estimated period
    before the first extremum belonging to the sustained oscillatory
    regime.

    Parameters
    ----------
    time : array_like
        Time coordinates.
    values : array_like
        Time-series values.
    prominence : float
        Minimum prominence used to detect maxima and minima.
    min_cycles : int, optional
        Minimum number of complete cycles required. Default is 3.
    period_tolerance : float, optional
        Maximum fractional deviation of individual full periods from
        the median period. Default is 0.25.
    min_extrema : int, optional
        Required number of alternating extrema, overriding min_cycles.
        Must be at least 4 so both maxima and minima define a full period.

    Returns
    -------
    t_onset : float
        Estimated onset time.
    period : float
        Estimated oscillation period.
    extrema : ndarray
        Indices of the extrema used to identify the oscillatory
        regime.
    """
    time = np.asarray(time)
    values = np.asarray(values)

    if time.ndim != 1 or values.ndim != 1:
        raise ValueError("time and values must be one-dimensional.")

    if len(time) != len(values):
        raise ValueError("time and values must have equal length.")

    if len(time) < 2 or np.any(np.diff(time) <= 0):
        raise ValueError("time must be strictly increasing.")

    # Find prominent maxima and minima.
    maxima, _ = find_peaks(values, prominence=prominence)
    minima, _ = find_peaks(-values, prominence=prominence)

    # Combine and sort all extrema.
    indices = np.concatenate((maxima, minima))
    types = np.concatenate(
        (
            np.ones(len(maxima), dtype=int),
            -np.ones(len(minima), dtype=int),
        )
    )

    order = np.argsort(indices)
    indices = indices[order]
    types = types[order]

    # Three cycles require seven alternating extrema:
    #
    # max min max min max min max
    #
    # or vice versa.
    n_required = 2 * min_cycles + 1
    if min_extrema is not None:
        if not isinstance(min_extrema, (int, np.integer)) or min_extrema < 4:
            raise ValueError("min_extrema must be an integer of at least 4.")
        n_required = min_extrema

    if len(indices) < n_required:
        raise ValueError(
            "Too few extrema to identify sustained oscillations."
        )

    for start in range(len(indices) - n_required + 1):
        stop = start + n_required

        candidate_indices = indices[start:stop]
        candidate_types = types[start:stop]

        # Maxima and minima must alternate.
        if not np.all(candidate_types[1:] == -candidate_types[:-1]):
            continue

        candidate_maxima = candidate_indices[
            candidate_types == 1
        ]
        candidate_minima = candidate_indices[
            candidate_types == -1
        ]

        # Calculate full periods, never half-periods.
        max_periods = np.diff(time[candidate_maxima])
        min_periods = np.diff(time[candidate_minima])

        periods = np.concatenate((max_periods, min_periods))

        if len(periods) == 0:
            continue

        period = np.median(periods)

        # Require all complete periods to be reasonably consistent.
        fractional_deviation = np.abs(periods - period) / period

        if np.any(fractional_deviation > period_tolerance):
            continue

        # First extremum belonging to the sustained oscillation.
        t_first_extremum = time[candidate_indices[0]]

        # Phase-based estimate of oscillation onset.
        t_onset = t_first_extremum - period / 2.0

        return t_onset, period, candidate_indices

    raise ValueError("No sustained oscillatory regime found.")

def first_zero(time, values):
    """Return the first positive-to-negative zero crossing."""

    crossings = np.where(
        (values[:-1] > 0) & (values[1:] <= 0)
    )[0]

    if len(crossings) == 0:
        raise ValueError("No positive-to-negative reversal found.")

    i = crossings[0]

    return (
        time[i]
        - values[i]
        * (time[i + 1] - time[i])
        / (values[i + 1] - values[i])
    )

def time_first_minimum(time, values, prominence=None):
    """Return the time of the first local minimum in the time series."""

    minima, _ = find_peaks(-values, prominence=prominence)

    if len(minima) == 0:
        raise ValueError("No reversal found.")

    return time[minima[0]]
