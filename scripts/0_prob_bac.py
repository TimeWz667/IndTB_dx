from sim.healthcare.system import *


if __name__ == '__main__':
    import json

    p0 = {
        'sens_ssm': 0.64,
        'spec_ssm': 0.98,
        'sens_xpert': 0.85,
        'sens_xpert_ss-': 0.64,
        'spec_xpert': 0.98,
        'p_ava_xpert_pub': 0.2,
        'p_ava_ssm_pub': 0.8,
        'p_ava_xpert_eng': 0.3,
        'p_loss_sputum': 0.15,
        'p_loss_swab': 0.02,
        'sens_cdx': 0.7,
        'spec_cdx': 0.95,
        'p_csi_pub': 0.483,
        'p_csi_ppm': 0.6,
        'dur_pri': 0.7,
        'dur_pub': 0.5,
        'p_txi_pub': 0.9,
        'p_txi_eng': 0.8,
        'p_txi_pri': 0.8,
        'p_refer_i2u': 0.3,
    }

    system = get_system(p0, has_cdx=False)

    tp_bac = list()
    test_tb_ssm = list()
    test_tb_xpert = list()
    mb_tb = list()

    fp_bac = list()
    test_nontb_ssm = list()
    test_nontb_xpert = list()
    mb_nontb = list()

    print('---------------------------------------------')
    for alg in system.Public.Algorithms:
        res0 = alg.dx(1, 0)
        res1 = alg.dx(0, 1)

        tp_bac.append(res0.TruePos)
        test_tb_ssm.append(res0['N_test_SSM'])
        test_tb_xpert.append((res0['N_test_Xpert_ss-'] + res0['N_test_Xpert']))
        mb_tb.append(res0['N_MisBacLTFU'])

        fp_bac.append(res1.FalsePos)
        test_nontb_ssm.append(res1['N_test_SSM'])
        test_nontb_xpert.append((res1['N_test_Xpert_ss-'] + res1['N_test_Xpert']))
        mb_nontb.append(res1['N_MisBacLTFU'])

    print(tp_bac)
    print(fp_bac)
    print(test_tb_ssm)
    print(test_nontb_ssm)
    print(mb_tb)
    print(mb_nontb)

    js = dict(
        tp_bac=tp_bac,
        fp_bac=fp_bac,
        test_tb_ssm=test_tb_ssm,
        test_tb_xpert=test_tb_xpert,
        test_nontb_ssm=test_nontb_ssm,
        test_nontb_xpert=test_nontb_xpert,
        mb_tb=mb_tb,
        mb_nontb=mb_nontb
    )

    with open('../data/pre_cdx.json', 'w') as f:
        json.dump(js, f)


    # tp_bac = c(0.73984, 0.544, 0.7225, 0),
    # fp_bac = c(0.0336, 0.017, 0.017, 0),
    # test_tb_ssm = c(0.85, 0.85, 0, 0),
    # test_tb_naat = c(0.306, 0, 0.85, 0),
    # test_nontb_ssm = c(0.85, 0.85, 0, 0),
    # test_nontb_naat = c(0.833, 0, 0.85, 0),

