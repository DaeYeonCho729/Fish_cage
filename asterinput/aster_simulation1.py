import sys
import os
import numpy as np

cwd = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cwd, "module"))

import morison as mori
import screen as scr  
import aster_module

DEBUT(
    PAR_LOT='NON',
    IGNORE_ALARM=("SUPERVIS_25", "DISCRETE_26", "UTILITAI8_56")
)

INCLUDE(UNITE=90)
INCLUDE(UNITE=91)


dt = materials["simulation"]["dt"]
duration = materials["simulation"]["duration"]

r = (meshinfo["Numerical_half_mesh_size"]/materials["Net"]["Mesh_length"])
rr = r **0.5
current_input = materials['simulation']['current']

itimes=int(duration/dt)
tend = itimes*dt

nodenumber = meshinfo["numberofNodes"]
l = ["None"] * (nodenumber + 1) 

S_a = meshinfo["S_a"]

netting = scr.ScreenModel(
    rr,
    meshinfo["surfs_netting"],     
    materials["Net"]["Sn"],        
    materials["Net"]["twine_dia"], 
    materials["Net"]["RHO"],      
    lines_netting = meshinfo["Lines_netting"],  
    wake_origin = meshinfo["wake_origin"] 
)

floater_morison = mori.MorisonModel(
    meshinfo["Lines_pipe_top"],
    materials["Floater"]["Floaterdiameter"],
    materials["Floater"]["RHO"],
    materials["Floater"]["Floaterthickness"] * 0.5
)

bottomring_morison = mori.MorisonModel(
    meshinfo["Lines_bottom_ring"],
    materials["Bottom_ring"]["Br_diameter"],
    materials["Bottom_ring"]["RHO"],
    materials["Bottom_ring"]["Br_diameter"] * 0.5
)

siderope_morison = mori.MorisonModel(
    meshinfo["Lines_side_rope"],
    materials["Side_ropes"]["SideRope_diameter"],
    materials["Side_ropes"]["RHO"],
    materials["Side_ropes"]["SideRope_diameter"] * 0.5
)

bridle_morison = mori.MorisonModel(
    meshinfo["Lines_bridle"],
    materials["Bridle"]["Bridle_diameter"],
    materials["Bridle"]["RHO"],
    materials["Bridle"]["Bridle_diameter"] * 0.5
)

buoyline_morison = mori.MorisonModel(
    meshinfo["Lines_buoy"],
    materials["Buoyline1"]["Buoy1_diameter"],
    materials["Buoyline1"]["RHO"],
    materials["Buoyline1"]["Buoy1_diameter"] * 0.5
)

anchora_morison = mori.MorisonModel(
    meshinfo["Lines_anchor_A"],
    materials["Anchor_A"]["AnchorA_diameter"],
    materials["Anchor_A"]["RHO"],
    materials["Anchor_A"]["AnchorA_diameter"] * 0.5
)

anchorb_morison = mori.MorisonModel(
    meshinfo["Lines_anchor_B"],
    materials["Anchor_B"]["AnchorB_diameter"],
    materials["Anchor_B"]["RHO"],
    materials["Anchor_B"]["AnchorB_diameter"] * 0.5
)

mooring_morison = mori.MorisonModel(
    meshinfo["Lines_mooring_frame"],
    materials["Mooring_frame"]["MoorFrame_diameter"],
    materials["Mooring_frame"]["RHO"],
    materials["Mooring_frame"]["MoorFrame_diameter"] * 0.5
)

braket_morison = mori.MorisonModel(
    meshinfo["Line_braket"],
    materials["Bracket"]["Brkt_diameter"],
    materials["Bracket"]["RHO"],
    materials["Bracket"]["Brkt_diameter"] * 0.5
)

n_net  = len(netting.output_hydro_element())
n_flo  = len(floater_morison.output_hydro_element())
n_btr  = len(bottomring_morison.output_hydro_element())
n_sir  = len(siderope_morison.output_hydro_element())
n_bri  = len(bridle_morison.output_hydro_element())
n_buyl = len(buoyline_morison.output_hydro_element())
n_aa   = len(anchora_morison.output_hydro_element())
n_ab   = len(anchorb_morison.output_hydro_element())
n_moor = len(mooring_morison.output_hydro_element())
n_bra  = len(braket_morison.output_hydro_element())

i0 = 0
s_net  = slice(i0, i0 + n_net );  i0 += n_net
s_flo  = slice(i0, i0 + n_flo );  i0 += n_flo
s_btr  = slice(i0, i0 + n_btr );  i0 += n_btr
s_sir  = slice(i0, i0 + n_sir );  i0 += n_sir
s_bri  = slice(i0, i0 + n_bri );  i0 += n_bri
s_buyl = slice(i0, i0 + n_buyl);  i0 += n_buyl
s_aa   = slice(i0, i0 + n_aa  );  i0 += n_aa
s_ab   = slice(i0, i0 + n_ab  );  i0 += n_ab
s_moor = slice(i0, i0 + n_moor);  i0 += n_moor
s_bra  = slice(i0, i0 + n_bra );  i0 += n_bra

NO_cage = i0

#v2024/Fish_cage/asteroutput/hydro_elements.txt
out_dir = os.path.join(cwd, "asteroutput")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "hydro_elements.txt")

with open(out_path, "w", encoding="utf-8") as file:
    # 1) netting
    file.write("netting\n")
    file.write(str(netting.output_hydro_element()))
    file.write("\n")

# read mesh
mesh = LIRE_MAILLAGE(
    FORMAT="IDEAS",   
    UNITE=20,
    INFO=1,
    VERI_MAIL=_F(VERIF='OUI', APLAT=0.001),
)

affe = []

if not S_a:
    affe.append(
        _F(
            GROUP_MA=('NETTWIN', 'SIDEROPE', 'BUOYLN1', 'ANCHL1', 'BRIDLE', 'ANCHL2', 'MOORFRM'),
            MODELISATION='CABLE',
            PHENOMENE='MECANIQUE'
        )
    )
    affe.append(
        _F(
            GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING'),
            MODELISATION='POU_D_E',
            PHENOMENE='MECANIQUE'
        )
    )
else:
    affe.append(
        _F(
            GROUP_MA=('NETTWIN', 'BOTRING', 'BUOYLN1', 'ANCHL1', 'BRIDLE', 'ANCHL2', 'MOORFRM'),
            MODELISATION='CABLE',
            PHENOMENE='MECANIQUE'
        )
    )
    affe.append(
        _F(
            GROUP_MA=('FLOATER', 'BRACKET'),
            MODELISATION='POU_D_E',
            PHENOMENE='MECANIQUE'
        )
    )

model = AFFE_MODELE(
    AFFE=tuple(affe),
    MAILLAGE=mesh
)

if not S_a:
    elemprop = AFFE_CARA_ELEM( 
        CABLE=(
            _F(
            GROUP_MA=('NETTWIN',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Net"]["twine_dia"] * rr , 2)
        ),
            _F(
            GROUP_MA=('BUOYLN1',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Buoyline1"]["Buoy1_diameter"], 2)
        ),
            _F(
            GROUP_MA=('ANCHL1',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Anchor_A"]["AnchorA_diameter"], 2)
        ),
            _F(
            GROUP_MA=('BRIDLE',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Bridle"]["Bridle_diameter"], 2)
        ),
            _F(
            GROUP_MA=('ANCHL2',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Anchor_B"]["AnchorB_diameter"], 2)
        ),
            _F(
            GROUP_MA=('MOORFRM',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Mooring_frame"]["MoorFrame_diameter"], 2)
        ),
            _F(
            GROUP_MA=('SIDEROPE',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Side_ropes"]["SideRope_diameter"], 2)
        ),
        ),
        POUTRE=(
            _F(
                GROUP_MA=('FLOATER', ),
                SECTION='CERCLE',
                CARA=('R', 'EP'),
                VALE=(
                    materials["Floater"]["Floaterdiameter"] / 2.0,
                    materials["Floater"]["Floaterthickness"]
                )
            ),
            _F(
                GROUP_MA=('BOTRING',),
                SECTION='CERCLE',
                CARA=('R',),
                VALE=(
                    materials["Bottom_ring"]["Br_diameter"] / 2.0,
                )
            ),
            _F(
                GROUP_MA=('BRACKET'),
                SECTION='CERCLE',
                CARA=('R',),
                VALE=(
                    materials["Bracket"]["Brkt_diameter"] / 2.0,
                )
            ),
        ),
        MODELE=model
    )
else:
    elemprop = AFFE_CARA_ELEM( 
        CABLE=(_F(
            GROUP_MA=('NETTWIN','BOTRING'),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Net"]["twine_dia"] * rr , 2)
        ),
            _F(
            GROUP_MA=('BUOYLN1',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Buoyline1"]["Buoy1_diameter"], 2)
        ),
            _F(
            GROUP_MA=('ANCHL1',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Anchor_A"]["AnchorA_diameter"], 2)
        ),
            _F(
            GROUP_MA=('BRIDLE',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Bridle"]["Bridle_diameter"], 2)
        ),
            _F(
            GROUP_MA=('ANCHL2',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Anchor_B"]["AnchorB_diameter"], 2)
        ),
            _F(
            GROUP_MA=('MOORFRM',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Mooring_frame"]["MoorFrame_diameter"], 2)
        ),
            _F(
            GROUP_MA=('SIDEROPE',),
            N_INIT=10.0,
            SECTION=0.25 * np.pi * pow(materials["Side_ropes"]["SideRope_diameter"], 2)
        ),
        ),
        POUTRE=(
            _F(
                GROUP_MA=('FLOATER', 'BRACKET'),
                SECTION='CERCLE',
                CARA=('R', 'EP'),
                VALE=(
                    materials["Floater"]["Floaterdiameter"] / 2.0,
                    materials["Floater"]["Floaterthickness"]
                )
            ),
        ),
        MODELE=model
    )      

net = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Net"]["NetYoungModule"],
        NU  = 0.3,
        RHO = materials["Net"]["RHO"]
    )
)

buoyline = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Buoyline1"]["Buoy1YoungModule"],
        NU  = 0.3,
        RHO = materials["Buoyline1"]["RHO"]
    )
)

anchorline1 = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Anchor_A"]["AnchorAYoungModule"],
        NU  = 0.3,
        RHO = materials["Anchor_A"]["RHO"]
    )
)

anchorline2 = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Anchor_B"]["AnchorBYoungModule"],
        NU  = 0.3,
        RHO = materials["Anchor_B"]["RHO"]
    )
)

bridle = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Bridle"]["BridleYoungModule"],
        NU  = 0.3,
        RHO = materials["Bridle"]["RHO"]
    )
)

mooringframe = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Mooring_frame"]["MoorFrameYoungModule"],
        NU  = 0.3,
        RHO = materials["Mooring_frame"]["RHO"]
    )
)

siderope = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Side_ropes"]["SideRopeYoungModule"],
        NU  = 0.3,
        RHO = materials["Side_ropes"]["RHO"]
    )
)

floater = DEFI_MATERIAU(
    ELAS=_F(
        
        E   = materials["Floater"]["FloaterYoungModule"],
        NU  = 0.3,
        RHO = materials["Floater"]["RHO"]
    )
)

braket = DEFI_MATERIAU(
    ELAS=_F(
        E   = materials["Bracket"]["BrktYoungModule"],
        NU  = 0.3,
        RHO = materials["Bracket"]["RHO"]
    )
)

if not S_a:
    bottomring = DEFI_MATERIAU(
        ELAS=_F(
            E   = materials["Bottom_ring"]["BrYoungModule"],
            NU  = 0.3,
            RHO = materials["Bottom_ring"]["RHO"]
        )
    )
else:
    bottomring = DEFI_MATERIAU(
        CABLE=_F(EC_SUR_E=0.0001),
        ELAS=_F(
            E   = materials["Net"]["NetYoungModule"],
            NU  = 0.3,
            RHO = materials["Net"]["RHO"]
        )
    )

if not S_a:
    fieldmat = AFFE_MATERIAU(
        AFFE=(
            _F(GROUP_MA=('NETTWIN',),  MATER=net),
            _F(GROUP_MA=('FLOATER',),  MATER=floater),
            _F(GROUP_MA=('BRACKET',),  MATER=braket),
            _F(GROUP_MA=('BOTRING',),  MATER=bottomring),
            _F(GROUP_MA=('MOORFRM',),  MATER=mooringframe),
            _F(GROUP_MA=('BRIDLE',),  MATER=bridle),
            _F(GROUP_MA=('SIDEROPE',),  MATER=siderope),
            _F(GROUP_MA=('BUOYLN1',),  MATER=buoyline),
            _F(GROUP_MA=('ANCHL1',),  MATER=anchorline1),
            _F(GROUP_MA=('ANCHL2',),  MATER=anchorline2),
        ),
        MODELE=model
    )

else:
    fieldmat = AFFE_MATERIAU(
        AFFE=(
            _F(GROUP_MA=('NETTWIN','BOTRING'),  MATER=net), 
            _F(GROUP_MA=('FLOATER',),           MATER=floater),
            _F(GROUP_MA=('BRACKET',),           MATER=braket),
            _F(GROUP_MA=('BOTRING',),  MATER=bottomring),
            _F(GROUP_MA=('MOORFRM',),  MATER=mooringframe),
            _F(GROUP_MA=('BRIDLE',),  MATER=bridle),
            _F(GROUP_MA=('SIDEROPE',),  MATER=siderope),
            _F(GROUP_MA=('BUOYLN1',),  MATER=buoyline),
            _F(GROUP_MA=('ANCHL1',),  MATER=anchorline1),
            _F(GROUP_MA=('ANCHL2',),  MATER=anchorline2),
        ),
        MODELE=model
    )

# if not S_a:
#     #gravity 
#     gF = AFFE_CHAR_MECA(
#         PESANTEUR=_F(
#             DIRECTION=(0.0, 0.0, -1.0),
#             GRAVITE=9.81,
#             GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING','MOORFRM','BRIDLE','SIDEROPE','BUOYLN1','ANCHL1','ANCHL2',),
#         ),
#         MODELE=model
#     )
# else:
#     #gravity
#     gF = AFFE_CHAR_MECA(
#         PESANTEUR=_F(
#             DIRECTION=(0.0, 0.0, -1.0),
#             GRAVITE=9.81,
#             GROUP_MA=('FLOATER', 'BRACKET','MOORFRM','BRIDLE','SIDEROPE','BUOYLN1','ANCHL1','ANCHL2',),
#         ),
#         MODELE=model
#     )

#floater fix
fixed = AFFE_CHAR_MECA(
    DDL_IMPO=_F(
        GROUP_NO=('ANCHORN',),   # floater node → FLONODE
        LIAISON='ENCASTRE'
    ),
    MODELE=model
)


if materials["Net"]["bottom_net_center"]:
    load_list1 = [
        _F(
            GROUP_NO=('BOTNCEN',),
            FX=0,
            FY=0,
            FZ=-materials["Net"]["BottomNet_cen_weight"],
        )
    ]
else:
    load_list1 = [
        _F(
            GROUP_NO=('BOTNCEN',),
            FX=0,
            FY=0,
            FZ=0,
        ),
    ]


sinkF1 = AFFE_CHAR_MECA(
    FORCE_NODALE=tuple(load_list1),
    MODELE=model
)

# if sinker for bottomring node -> weight force 
if S_a:
    load_list2 = [
        _F(
            GROUP_NO=('BRNODE',),
            FX=0,
            FY=0,
            FZ=-materials["Sinker"]["single_sinker_weight"],
        )
    ]
else:
    load_list2 = [
        _F(
            GROUP_NO=('BRNODE',),
            FX=0,
            FY=0,
            FZ=0,
        ),
    ]

sinkF2 = AFFE_CHAR_MECA(
    FORCE_NODALE=tuple(load_list2),
    MODELE=model
)

# if sinker for bottomring node -> weight force 

load_list3 = [
    _F(
        GROUP_NO=('BUOY1',),
        FX=0,
        FY=0,
        FZ=400,
    )
]


sinkF3 = AFFE_CHAR_MECA(
    FORCE_NODALE=tuple(load_list3),
    MODELE=model
)


load_list4 = [
    _F(
        GROUP_NO=('ANCHORN',),
        FX=0,
        FY=0,
        FZ=-3000,
    )
]


sinkF4 = AFFE_CHAR_MECA(
    FORCE_NODALE=tuple(load_list3),
    MODELE=model
)

listr = DEFI_LIST_REEL(DEBUT=0.0,
                       INTERVALLE=_F(JUSQU_A=tend,PAS=dt))

times = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=listr,PAS_MINI=1e-8),
                       METHODE='AUTO')

for k in range(0,itimes):
    INCLUDE(UNITE=93,INFO=0) #UNITE=93 : aster_simulation2

# stat1 = CALC_CHAMP(
#     RESULTAT=resn,
#     GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING'),
#     CONTRAINTE=('EFGE_ELGA', 'EFGE_ELNO', 'SIPO_ELNO'),
# )

# stat2 = CALC_CHAMP(
#     RESULTAT=resn,
#     GROUP_MA=('NETTWIN',),
#     CONTRAINTE=('EFGE_ELGA', 'EFGE_ELNO'),
# )

#  #reaction force 
# stat3 = CALC_CHAMP(RESULTAT=resn,
#                   CONTRAINTE=('SIEF_ELNO',
#                               ),
#                   FORCE=('REAC_NODA', ),
#                   )

# IMPR_RESU(
#     FORMAT='MED',
#     RESU=(
#         _F(
#             RESULTAT=resn,
#             GROUP_MA=('NETTWIN','FLOATER','BRACKET','BOTRING'),
#             NOM_CHAM=('DEPL','SIEF_ELGA'),
#             LIST_INST=listr,
#             TOUT_CMP='OUI',
#         ),
#         _F(
#             RESULTAT=stat1,
#             GROUP_MA=('FLOATER','BRACKET','BOTRING'),
#             NOM_CHAM=('EFGE_ELGA','EFGE_ELNO','SIPO_ELNO'),
#             LIST_INST=listr,
#         ),
#         _F(
#             RESULTAT=stat2,
#             GROUP_MA=('NETTWIN',),
#             NOM_CHAM=('EFGE_ELGA','EFGE_ELNO'),
#             LIST_INST=listr,
#         ),
#     ),
#     UNITE=80
# )


# reac1 = POST_RELEVE_T(
#     ACTION=_F(
#         GROUP_NO=('FLONODE',),  # floater node  FLONODE
#         INTITULE='sum reactions',
#         MOMENT=('DRX', 'DRY', 'DRZ'),
#         NOM_CHAM=('REAC_NODA',),
#         OPERATION=('EXTRACTION',),
#         POINT=(0.0, 0.0, 0.0),
#         RESULTANTE=('DX', 'DY', 'DZ'),
#         LIST_INST=listr,     
#         RESULTAT=stat3,
#     )
# )                         


# IMPR_TABLE(FORMAT_R='1PE12.3',
#            TABLE=reac1,
#            UNITE=10)                    


FIN()