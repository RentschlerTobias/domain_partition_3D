"""Stand-in plugins for the registry and contract tests."""

from hexdp.contracts import ParamSpec, Score


class DummyStrategy:
    name = "dummy"
    params = {"tol": ParamSpec(0.05, "tolerance", lower=0)}

    def __init__(self, tol: float = 0.05):
        self.tol = tol

    def applicable(self, geo, features):
        return Score(0.1, ["dummy"])

    def run(self, geo, features, params, ctx):
        raise NotImplementedError


class NoParams:
    name = "noparams"

    def applicable(self, geo, features):
        return Score.no("never")
