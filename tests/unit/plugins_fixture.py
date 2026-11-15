"""Stand-in plugin for the registry tests."""


class DummyStrategy:
    name = "dummy"

    def __init__(self, tol: float = 0.05):
        self.tol = tol
