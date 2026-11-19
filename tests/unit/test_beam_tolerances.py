"""The beam's two Hausdorff tolerances: relative (fraction of the bbox diagonal) is the default."""
from hexdp.ops.beam_collapse import DEFAULT_ABSOLUTE, DEFAULT_RELATIVE, tolerances


def test_defaults_are_relative():
    tol, raw = tolerances(None, None, 100.0, True)
    assert (tol, raw) == (DEFAULT_RELATIVE[0] * 100.0, DEFAULT_RELATIVE[1] * 100.0)


def test_defaults_in_absolute_mode_keep_the_old_numbers():
    assert tolerances(None, None, 100.0, False) == DEFAULT_ABSOLUTE


def test_given_values_are_scaled_only_in_relative_mode():
    assert tolerances(0.01, 0.2, 50.0, True) == (0.5, 10.0)
    assert tolerances(0.01, 0.2, 50.0, False) == (0.01, 0.2)


def test_the_cli_defaults_to_relative_and_has_an_absolute_switch():
    """The parser lives in main(), so read it from the source: --relative is the default, --absolute turns it off."""
    import inspect

    from hexdp.ops import beam_collapse

    src = inspect.getsource(beam_collapse)
    assert 'dest="relative", action="store_true", default=True' in src
    assert 'dest="relative", action="store_false"' in src
