element = [
'Fe',
'Co',
'Nb',
'Tc',
'V',
'Zr',
'Ag',
'Cu',
'Pd',
'Rh',
'Ru',
'Sc',
'Ti',
'Y',
'Ni',
'Cr',
'Mo',
'Mn',
'Au',
'Hf',
'Ir',
'Os',
'Pt',
'Re',
'Ta',
'W',
'Cu',
]

fe = [
5.158732468,
4.48889952,
1.49547267,
4.575253135,
2.48340383,
-0.50558209,
7.27840966,
None,
6.08428732,
6.0129126,
5.454666625,
-1.71209207,
0,
-1.377292445,
5.58351351,
None,
3.26196065,
None,
8.18370044,
-0.915791655,
7.06715516,
6.307231145,
6.8575222,
5.055205595,
1.68971541,
3.99039433,
6.55473439999332,
]

plus_4_fe = {
        'Sc':0.630560999999,
'Ti':-4.999990324e-07,
'V':2.8244667,
'Fe':8.4013777,
'Co':7.4183352,
'Ni':9.2910159,
'Cu':10.157252,
'Y':1.1001678,
'Zr':-0.748096899998,
'Nb':1.3983322,
'Mo':3.9490165,
'Tc':5.877902,
'Ru':7.5010117,
'Rh':8.7066408,
'Pd':9.8778824,
'Ag':11.0927297,
'Hf':-1.2150706,
'Ta':1.3929598,
'W':4.3553636,
'Re':5.8966863,
'Os':7.6889562,
'Ir':9.1189437,
'Pt':10.084685,
'Au':11.2244594,
}



column = [
8,
9,
5,
7,
5,
4,
11,
11,
10,
9,
8,
3,
4,
3,
10,
6,
6,
7,
11,
4,
9,
8,
10,
7,
5,
6,
]


fe_dict = {a:b for a, b in zip(element, fe)}
#with open('../data/data/formation_energy.csv', 'w') as f:
#    for key, value in fe_dict.items():
#        if value is None:
#            continue
#        f.write(key + ',' + str(value) + '\n')

column_dict = {a:b for a, b in zip(element, column)}


cohesive_energies = {
'Fe':4.28,
'Co':4.39,
'Nb':7.57,
'Tc':6.85,
'V':5.31,
'Zr':6.25,
'Ag':2.95,
'Cu':3.49,
'Pd':3.89,
'Rh':5.75,
'Ru':6.74,
'Sc':3.90,
'Ti':4.85,
'Y':4.37,
'Ni':4.44,
'Cr':4.10,
'Mo':6.82,
'Mn':2.92,
'Au':3.81,
'Hf':6.44,
'Ir':6.94,
'Os':8.17,
'Pt':5.84,
'Re':8.03,
'Ta':8.10,
'W':8.90,

}


d_cohesive = {
'Fe':476,
'Co':315 ,
'Nb':1156 ,
'Tc':1060,
'V':800,
'Zr':921,
'Ag':0,
'Cu':0,
'Pd':144,
'Rh':460,
'Ru':807,
'Sc':468,
'Ti':659,
'Y':612,
'Ni':136,
'Cr':788,
'Mo':1197,
'Mn':637,
'Au':0,
'Hf':1023,
'Ir':676,
'Os':1061,
'Pt':287,
'Re':1278,
'Ta':1303,
'W':1378,
}

s_cohesive = {
'Fe':422,
'Co':424 ,
'Nb':222,
'Tc':328,
'V':306,
'Zr':161,
'Ag':249,
'Cu':422,
'Pd':235,
'Rh':282,
'Ru':335,
'Sc':154,
'Ti':228,
'Y':117,
'Ni':398,
'Cr':359,
'Mo':289,
'Mn':406,
'Au':452,
'Hf':291,
'Ir':587,
'Os':652,
'Pt':508,
'Re':604,
'Ta':407,
'W':527,

}


electronegativity = {
'Fe':1.83,
'Co':1.88,
'Nb':1.60,
'Tc':1.90,
'V':1.63,
'Zr':1.33,
'Ag':1.93,
'Cu':1.90,
'Pd':2.20,
'Rh':2.28,
'Ru':2.20,
'Sc':1.36,
'Ti':1.54,
'Y':1.22,
'Ni':1.91,
'Cr':1.66,
'Mo':2.16,
'Mn':1.55,
'Au':2.54,
'Hf':1.30,
'Ir':2.20,
'Os':2.2,
'Pt':2.28,
'Re':1.9,
'Ta':1.5,
'W':2.36,

}

N2_engs = [
-0.900475549,
-0.275275798,
-0.573250298,
-1.028873298,
-0.509580198,
-0.437712498,
-0.152380398,
-0.207399798,
-0.149891798,
-0.479393598,
-1.247303798,
-0.355808998,
-0.335471498,
-0.367079998,
-0.215497298,
None,
-0.649177998,
None,
-0.128035098,
-0.443044898,
-0.987119298,
-1.583956098,
-0.121257398,
-1.188556598,
-0.546712598,
-0.677330198,
]

N2_engs_dict = {a:b for a, b in zip(element, N2_engs)}

N2H_engs = [
0.9627901614,
0.0920799095,
-0.5031168905,
-0.7638951905,
-0.3414724905,
-0.0476147905,
1.75579361,
1.49092071,
1.209224809,
-0.2206814905,
-0.1681130905,
0.6845924095,
0.1167635095,
-0.1192273905,
0.7911753095,
None,
-0.7515885905,
None,
1.644681909,
-0.1640615905,
-0.6723059905,
-0.9212194905,
1.03340891,
-1.124197491,
-0.5697134905,
-0.9835309905,
]

N2H_engs_dict = {a:b for a, b in zip(element, N2H_engs)}

NH2_engs_dict = {'Sc': -0.4282966984023262, 'Ti': -0.8185341983984729, 'V': -0.8612199983976354, 'Co': -0.5669012984066542, 'Ni': 0.11657690159388201, 'Cu': 0.9528160015984866, 'Y': 0.023393401600125596, 'Zr': -1.0885280984070307, 'Nb': -1.5198739983965681, 'Mo': -1.0271336984073316, 'Tc': -0.6943065984009338, 'Ru': -0.7558839984043773, 'Rh': -0.8842827984039361, 'Pd': 0.43984160159795627, 'Ag': 1.2892380015915221, 'Hf': -1.2616272984009158, 'Ta': -1.686536898406692, 'W': -1.426843298402153, 'Re': -1.0328282984041182, 'Os': -1.2625424984066087, 'Ir': -1.3418939984034552, 'Pt': 0.07911950159957248, 'Au': 0.9211964015926895}

d_band = {
'Fe':-0.92,
'Co':-1.17,
'Nb':1.41,
'Tc':-0.60,
'V':1.06,
'Zr':1.95,
'Ag':-4.30,
'Cu':-2.67,
'Pd':-1.83,
'Rh':-1.73,
'Ru':-1.41,
'Sc':None,
'Ti':1.50,
'Y':None,
'Ni':-1.29,
'Cr':0.16,
'Mo':0.35,
'Mn':0.07,
'Au':-3.56,
'Hf':2.47,
'Ir':-2.11,
'Os':None,
'Pt':-2.25,
'Re':-0.51,
'Ta':2.00,
'W':0.77,
}


plus_4_N2 = {
'Sc':-0.30651589999934004,
'Ti':-0.23107109999864406,
'V':-0.19632009999941147,
'Fe':-0.5072335999971074,
'Co':-0.1916634000040176,
'Ni':-0.18482649999668865,
'Cu':-0.1642661000037151,
'Y':-0.3506748999966476,
'Zr':-0.31588910000164105,
'Nb':-0.27760470000384885,
'Mo':-0.24728819999609186,
'Tc':-0.2699792000016714,
'Ru':-0.3903132999944319,
'Rh':-0.5318491000030008,
'Pd':-0.3949569999999767,
'Ag':-0.14717200000040975,
'Hf':-0.3259544999996251,
'Ta':-0.2844662000020435,
'W':-0.2541505999938636,
'Re':-0.26596700000254714,
'Os':-0.6112090999990869,
'Ir':-0.7796040999996876,
'Pt':-0.7219772999992529,
'Au':-0.11407799999460622,
}


plus_4_N2H = {
'Sc':1.4849088099976022,
'Ti':1.6917344099960125,
'V':1.9550116100030088,
'Ni':0.892486309999537,
'Y':1.7504768100001336,
'Zr':1.6045559099978952,
'Nb':1.5943271099998917,
'Rh':1.3678823100011748,
'Pd':0.8480230099988937,
'Hf':1.578166009999574,
'Ta':1.5124164100010482,
'W':1.5204165100034395,
'Os':1.1108214099980493,
'Ir':1.2430196099989068,
'Pt':0.9133459100006682,
'Au':1.3758149099989048,
}

bulk_fe_big = {"Fe": 8.017669599996225, "Co": 7.059392799992111, "Ir": 8.07600059999504, "Cu": 11.014249599992581, "Rh": 7.956369099995754, "Hf": -1.3314129000034427, "Mo": 3.5251062999964233, "Pd": 9.16208619999361, "Nb": 1.1514100999929724, "Au": 12.15569689999802, "Ru": 6.83728269999574, "V": 2.8201193999925636, "Ti": -2.5011104298755527e-12, "Re": 5.002647899993008, "Ni": 8.974257999991096, "Zr": -0.7819189000033475, "Os": 6.632604799997807, "Ag": 12.268288299995675, "Tc": 5.3298540999981014, "Sc": 0.75680099999704, "W": 3.6998673999955827, "Y": 1.436312099992847, "Pt": 9.071770799997466, "Ta": 0.9989691999976458}

