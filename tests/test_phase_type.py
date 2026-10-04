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
