import sys
import os
import numpy as np

cwd = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(cwd, "module"))

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

current_input = materials['simulation']['current']

itimes=int(duration/dt)
tend = itimes*dt

nodenumber = meshinfo["numberofNodes"]
l = ["None"] * (nodenumber + 1) 

netting = scr.ScreenModel(
    meshinfo["surfs_netting"],     
    materials["Net"]["Sn"],        
    materials["Net"]["twine_dia"], 
    materials["Net"]["RHO"],      
    lines_netting = meshinfo["Lines_netting"],  
    wake_origin = meshinfo["wake_origin"] 
)

NO_net = len(netting.output_hydro_element())

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

# define element type
model = AFFE_MODELE(
    AFFE=(
        _F(
            GROUP_MA=('NETTWIN',),        # net twine → NETTWIN
            MODELISATION='CABLE',
            PHENOMENE='MECANIQUE'
        ),
        _F(
            GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING'),  # floater, bracket, bottom ring
            MODELISATION='POU_D_E',
            PHENOMENE='MECANIQUE'
        ),
    ),
    MAILLAGE=mesh
)

elemprop = AFFE_CARA_ELEM(
    CABLE=_F(
        GROUP_MA=('NETTWIN',),
        N_INIT=10.0,
        SECTION=0.25 * np.pi * pow(materials["Net"]["twine_dia"] * 12 , 2) # need dwe 
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
        _F(
            GROUP_MA=('BOTRING',),
            SECTION='CERCLE',
            CARA=('R',),
            VALE=(
                materials["Bottom_ring"]["Br_diameter"] / 2.0,
            )
        ),
    ),
    MODELE=model
)

#define materials
net = DEFI_MATERIAU(
    CABLE=_F(EC_SUR_E=0.0001),
    ELAS=_F(
        E   = materials["Net"]["NetYoungModule"],
        NU  = 0.3,
        RHO = materials["Net"]["RHO"]
    )
)

floater = DEFI_MATERIAU(
    ELAS=_F(E=materials["Floater"]["FloaterYoungModule"],
    NU=0.3, 
    RHO=materials["Floater"]["RHO"]))

bottomring = DEFI_MATERIAU(
    ELAS=_F(E=materials["Bottom_ring"]["BrYoungModule"],
    NU=0.3, 
    RHO=materials["Bottom_ring"]["RHO"]))


fieldmat = AFFE_MATERIAU(
    AFFE=(
        _F(GROUP_MA=('NETTWIN',),  MATER=net),
        _F(GROUP_MA=('FLOATER',),  MATER=floater),
        _F(GROUP_MA=('BRACKET',),  MATER=floater),
        _F(GROUP_MA=('BOTRING',),  MATER=bottomring),
    ),
    MODELE=model
)


#add force
#gravity (bottomcenter)
gF = AFFE_CHAR_MECA(
    PESANTEUR=_F(
        DIRECTION=(0.0, 0.0, -1.0),
        GRAVITE=9.81,
        GROUP_MA=('NETTWIN', 'FLOATER', 'BRACKET', 'BOTRING'),
    ),
    MODELE=model
)

#floater fix
fixed = AFFE_CHAR_MECA(
    DDL_IMPO=_F(
        GROUP_NO=('FLONODE',),   # floater node → FLONODE
        LIAISON='ENCASTRE'
    ),
    MODELE=model
)


if materials["Net"]["bottom_net_center"]:
    load_list = [
        _F(
            GROUP_NO=('BOTNCEN',),
            FX=0,
            FY=0,
            FZ=-materials["Net"]["BottomNet_cen_weight"],
        )
    ]
else:
    load_list = [
        _F(
            GROUP_NO=('BOTNCEN',),
            FX=0,
            FY=0,
            FZ=0,
        ),
    ]


sinkF1 = AFFE_CHAR_MECA(
    FORCE_NODALE=tuple(load_list),
    MODELE=model
)


listr = DEFI_LIST_REEL(DEBUT=0.0,
                       INTERVALLE=_F(JUSQU_A=tend,PAS=dt))

times = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=listr,PAS_MINI=1e-8),
                       METHODE='AUTO')

for k in range(0,itimes):
    INCLUDE(UNITE=93,INFO=0) #UNITE=93 : aster_simulation2

stat1 = CALC_CHAMP(
    RESULTAT=resn,
    GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING'),
    CONTRAINTE=('EFGE_ELGA', 'EFGE_ELNO', 'SIPO_ELNO'),
)

stat2 = CALC_CHAMP(
    RESULTAT=resn,
    GROUP_MA=('NETTWIN',),
    CONTRAINTE=('EFGE_ELGA', 'EFGE_ELNO'),
)

#  #reaction force 
# stat3 = CALC_CHAMP(RESULTAT=resn,
#                   CONTRAINTE=('SIEF_ELNO',
#                               ),
#                   FORCE=('REAC_NODA', ),
#                   )

IMPR_RESU(
    FORMAT='MED',
    RESU=(
        _F(
            RESULTAT=resn,
            GROUP_MA=('NETTWIN','FLOATER','BRACKET','BOTRING'),
            NOM_CHAM=('DEPL','SIEF_ELGA'),
            LIST_INST=listr,
            TOUT_CMP='OUI',
        ),
        _F(
            RESULTAT=stat1,
            GROUP_MA=('FLOATER','BRACKET','BOTRING'),
            NOM_CHAM=('EFGE_ELGA','EFGE_ELNO','SIPO_ELNO'),
            LIST_INST=listr,
        ),
        _F(
            RESULTAT=stat2,
            GROUP_MA=('NETTWIN',),
            NOM_CHAM=('EFGE_ELGA','EFGE_ELNO'),
            LIST_INST=listr,
        ),
    ),
    UNITE=80
)


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