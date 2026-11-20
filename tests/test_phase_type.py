import numpy as np
import pytest

from hetload import PhaseType, marginal_jump_law, normalise_loadings


@pytest.fixture
def coxian():
    return PhaseType(alpha=np.array([0.7, 0.3]), T=np.array([[-3.0, 1.0], [0.0, -2.0]]))


def test_exponential():
    exp = PhaseType.exponential(2.0)
    assert exp.laplace(1.0) == pytest.approx(2.0 / 3.0)
    assert exp.mean() == pytest.approx(0.5)


@pytest.mark.parametrize("c", [0.25, 1.0, 3.0])
def test_scaling_matches_laplace(coxian, c):
    scaled = coxian.scale(c)
    for theta in [0.1, 1.0, 2.5 + 1.0j]:
        assert scaled.laplace(theta) == pytest.approx(coxian.laplace(c * theta))
    assert scaled.mean() == pytest.approx(c * coxian.mean())


def test_sample_mean(coxian):
    draws = coxian.sample(20_000, np.random.default_rng(0))
    assert draws.mean() == pytest.approx(coxian.mean(), rel=0.03)


@pytest.fixture(params=["exponential", "diagonal", "coxian", "cyclic"])
def law(request, coxian):
    if request.param == "exponential":
        return PhaseType.exponential(2.0)
    if request.param == "diagonal":
        return PhaseType.hyperexponential([0.4, 0.6], [1.0, 4.0])
    if request.param == "cyclic":
        return PhaseType(np.array([0.4, 0.6]), np.array([[-3.0, 2.0], [1.0, -2.0]]))
    return coxian


def test_sample_moments(law):
    first = np.linalg.solve(-law.T, np.ones(law.T.shape[0]))
    mean = law.alpha @ first
    variance = 2 * law.alpha @ np.linalg.solve(-law.T, first) - mean**2
    draws = law.sample(200_000, np.random.default_rng(42))
    assert draws.mean() == pytest.approx(mean, rel=0.015)
    assert draws.var() == pytest.approx(variance, rel=0.03)
    assert law.mean() == pytest.approx(mean)
    assert isinstance(law.mean(), float)


@pytest.mark.parametrize("shape", [(), (6,), (2, 3), (0,), (2, 0)])
def test_batched_laplace_matches_dense_solve(law, shape):
    theta = np.linspace(0.0, 3.0, int(np.prod(shape))).reshape(shape) + 0.3j
    expected = np.array(
        [
            law.alpha @ np.linalg.solve(t * np.eye(law.T.shape[0]) - law.T, law.exit_rates)
            for t in theta.flat
        ],
        dtype=complex,
    ).reshape(shape)
    actual = law.laplace(theta)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)
    if not shape:
        assert isinstance(actual, complex)
    else:
        assert actual.shape == shape
        assert np.iscomplexobj(actual)


def test_scalar_laplace_type(law):
    assert isinstance(law.laplace(1.0), complex)


def test_empty_sample(law):
    assert law.sample(0, np.random.default_rng(42)).shape == (0,)


def test_negative_sample_size(law):
    with pytest.raises(ValueError):
        law.sample(-1, np.random.default_rng(42))


def test_laplace_singular_matrix(coxian):
    with pytest.raises(np.linalg.LinAlgError):
        coxian.laplace(np.array([0.0, -3.0]))


def test_mutated_arrays_are_used():
    law = PhaseType.hyperexponential([0.4, 0.6], [1.0, 4.0])
    law.sample(10, np.random.default_rng(42))
    law.mean()
    law.laplace(1.0)
    law.T[0, 1] = 0.5
    law.alpha[:] = [1.0, 0.0]
    reference = PhaseType(law.alpha.copy(), law.T.copy())
    assert law.mean() == reference.mean()
    assert law.laplace(1.0) == reference.laplace(1.0)
    np.testing.assert_array_equal(
        law.sample(100, np.random.default_rng(42)),
        reference.sample(100, np.random.default_rng(42)),
    )


def test_marginal_law_homogeneous_matches_eq11():
    lam_j, eta, lam_c, gamma = 0.8, 2.0, 0.15, 3.5
    total, law = marginal_jump_law(lam_j, eta, lam_c, gamma)
    assert total == pytest.approx(lam_j + lam_c)
    theta = 0.7
    expected = (lam_j * eta / (eta + theta) + lam_c / (1.0 + gamma * theta)) / (lam_j + lam_c)
    assert law.laplace(theta) == pytest.approx(expected)


@pytest.mark.parametrize("omega", [0.5, 2.0])
def test_loading_scales_common_component(omega):
    lam_j, eta, lam_c, gamma = 0.8, 2.0, 0.15, 3.5
    _, loaded = marginal_jump_law(lam_j, eta, lam_c, gamma, omega)
    _, rescaled = marginal_jump_law(lam_j, eta, lam_c, omega * gamma)
    assert loaded.laplace(1.3) == pytest.approx(rescaled.laplace(1.3))


def test_normalise_loadings():
    omega = normalise_loadings([1.0, 2.0, 3.0])
    assert omega.mean() == pytest.approx(1.0)
    with pytest.raises(ValueError):
        normalise_loadings([1.0, -1.0])


@pytest.mark.parametrize(
    "law",
    [
        PhaseType.hyperexponential([0.3, 0.7], [2.0, 5.0]),
        PhaseType(np.array([1.0, 0.0]), np.array([[-2.0, 2.0], [0.0, -5.0]])),
    ],
)
def test_laplace_raises_at_pole(law):
    with pytest.raises(np.linalg.LinAlgError):
        law.laplace(-2.0)
    with pytest.raises(np.linalg.LinAlgError):
        law.laplace(np.array([1.0, -2.0]))
