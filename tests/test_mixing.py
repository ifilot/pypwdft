import numpy as np
import pytest

from pypwdft.mixing import DensityMixer


def test_pulay_accelerates_a_linear_fixed_point():
    target = np.array([1.0, 2.0, 4.0])
    contraction = np.array([0.85, 0.55, 0.20])

    linear = DensityMixer(np, method="linear", fraction=0.5, history=6)
    pulay = DensityMixer(np, method="pulay", fraction=0.5, history=6)
    linear_density = np.zeros(3)
    pulay_density = np.zeros(3)
    used_pulay = False
    for _ in range(6):
        linear_output = target + contraction * (linear_density - target)
        pulay_output = target + contraction * (pulay_density - target)
        linear_density = linear.update(linear_density, linear_output)
        pulay_density = pulay.update(pulay_density, pulay_output)
        used_pulay = used_pulay or pulay.last_step == "pulay"

    assert np.linalg.norm(pulay_density - target) < 0.6 * np.linalg.norm(
        linear_density - target
    )
    assert used_pulay


def test_mixer_preserves_the_requested_electron_density_mean():
    mixer = DensityMixer(
        np, method="pulay", fraction=0.5, history=4, target_mean=0.2
    )
    input_density = np.full((2, 2, 2), 0.2)
    output_density = np.linspace(0.1, 0.3, 8).reshape((2, 2, 2))

    mixed = mixer.update(input_density, output_density)

    assert np.mean(mixed) == pytest.approx(0.2)
    assert np.all(mixed >= 0)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"method": "broyden"}, "mixing must"),
        ({"fraction": 0}, "mixing_fraction"),
        ({"history": 1}, "mixing_history"),
    ],
)
def test_invalid_mixer_configuration(kwargs, message):
    with pytest.raises(ValueError, match=message):
        DensityMixer(np, **kwargs)
