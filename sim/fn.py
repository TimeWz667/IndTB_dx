from sim.ebm.dy import ModelBaseline
from sim.ebm.dy_risk import ModelRisk
from sim.inputs import load_inputs
from sim.ebm.intv import Intervention

__author__ = 'Chu-Chang Ku'
__all__ = ['load_model_baseline', 'load_model_powerlaw']


def load_model_baseline(root):
    inp = load_inputs(root)
    inp.Demography.HasMigration = False
    inp.Demography.set_year0(2000)
    return ModelBaseline(inp)


def load_model_powerlaw(root):
    inp = load_inputs(root)
    inp.Demography.HasMigration = False
    inp.Demography.set_year0(2000)
    return ModelRisk(inp)


if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    import matplotlib.pylab as plt

    model = load_model_powerlaw('../data')
    print(model)

    post = pd.read_csv('../pars/powerlaw/Post.csv').iloc[:, 1:]
    post = [dict(row) for _, row in post.iterrows()]
    post = [model.Inputs.Cascade.prepare_pars(p) for p in post]
    print(post[0])

    ys, ms, msg = model.simulate_to_fit(post[0], np.linspace(2016, 2020, 5))
    print(ms)

    y0 = ys.y[:, -1]

    intv1 = Intervention.parse_obj({
        'MassACF': {
            'Target': '20%',
            'Coverage': 0.2,
            'Sens_L': 0,
        }
    })

    intv2 = Intervention.parse_obj({
        'MassACF': {
            'Target': '20%',
            'Coverage': 0.2,
        }
    })

    _, ms0, _ = model.simulate_onward(y0, post[0])
    _, ms1, _ = model.simulate_onward(y0, post[0], intv1)
    _, ms2, _ = model.simulate_onward(y0, post[0], intv2)

    fig, axes = plt.subplots(4)

    ms0.IncR.plot(ax=axes[0])
    ms1.IncR.plot(ax=axes[0])
    ms2.IncR.plot(ax=axes[0])
    axes[0].set_title('Incidence')

    ms0.MorR.plot(ax=axes[1])
    ms1.MorR.plot(ax=axes[1])
    ms2.MorR.plot(ax=axes[1])
    axes[1].set_title('Mortality')

    ms0.YieldATBR.plot(ax=axes[2])
    ms1.YieldATBR.plot(ax=axes[2])
    ms2.YieldATBR.plot(ax=axes[2])
    axes[2].set_title('Yields, ATB')

    ms0.YieldOTR.plot(ax=axes[3])
    ms1.YieldOTR.plot(ax=axes[3])
    ms2.YieldOTR.plot(ax=axes[3])
    axes[3].set_title('Yields, OT')

    plt.show()
