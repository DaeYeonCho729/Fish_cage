import time
t0 = time.time()

save_interval = 0.01
save_step = int(save_interval / dt)

loadr = []
# loadr.append(_F(CHARGE=gF), )
loadr.append(_F(CHARGE=fixed), )
loadr.append(_F(CHARGE=sinkF1), )
loadr.append(_F(CHARGE=sinkF2), )
loadr.append(_F(CHARGE=sinkF3), )

if not S_a:
    comp = (
        _F(DEFORMATION='GROT_GDEP', GROUP_MA=('NETTWIN','MOORFRM','BRIDLE','SIDEROPE','BUOYLN1','ANCHL1','ANCHL2'), RELATION='CABLE'),
        _F(DEFORMATION='PETIT',     GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING'), RELATION='ELAS'),
    )
    obs_groups = ('NETTWIN','FLOATER','BRACKET','BOTRING','MOORFRM','BRIDLE','SIDEROPE','BUOYLN1','ANCHL1','ANCHL2',)
else:
    comp = (
        _F(DEFORMATION='GROT_GDEP', GROUP_MA=('NETTWIN','BOTRING','MOORFRM','BRIDLE','SIDEROPE','BUOYLN1','ANCHL1','ANCHL2'), RELATION='CABLE'),
        _F(DEFORMATION='PETIT',     GROUP_MA=('FLOATER','BRACKET'), RELATION='ELAS'),
    )
    obs_groups = ('NETTWIN','FLOATER','BRACKET','BOTRING','MOORFRM','BRIDLE','SIDEROPE','BUOYLN1','ANCHL1','ANCHL2',)



t1 = time.time()
if k == 0:
    resn = DYNA_NON_LINE(
        CARA_ELEM=elemprop,
        CHAM_MATER=fieldmat,
        COMPORTEMENT=comp,  
        CONVERGENCE=_F(ITER_GLOB_MAXI=100,
                       RESI_GLOB_RELA=1e-3),
        EXCIT=(loadr),
        OBSERVATION=_F(
            GROUP_MA=obs_groups,
            NOM_CHAM='DEPL',
            NOM_CMP=('DX', 'DY', 'DZ'),
            INST=k+dt,
            OBSE_ETAT_INIT='NON'
        ),
        SCHEMA_TEMPS=_F(FORMULATION='DEPLACEMENT',
                        SCHEMA="HHT",
                        ALPHA=-0.3),
        INCREMENT=_F(LIST_INST=times, INST_FIN=(1+k)*dt),
        MODELE=model)
    
else:
    Fnh = tuple(Fnh)
    force_list = []
    for i in range(nodenumber):
        force_list.append(_F(
            GROUP_NO=(f"N{i+1}",),   
            FX=Fnh[i][0],
            FY=Fnh[i][1],
            FZ=Fnh[i][2]
        ))

    load_dynamic = AFFE_CHAR_MECA(
        MODELE=model,
        FORCE_NODALE=tuple(force_list)
    )
    
    if len(loadr) == 4:
        loadr.append(_F(CHARGE=load_dynamic))
    else:
        loadr[4] = _F(CHARGE=load_dynamic)

    resn = DYNA_NON_LINE(CARA_ELEM=elemprop,
                         CHAM_MATER=fieldmat,
                         ETAT_INIT=_F(EVOL_NOLI=resn),
                         COMPORTEMENT=comp,
                         CONVERGENCE=_F(ITER_GLOB_MAXI=100,
                                        RESI_GLOB_RELA=1e-3),
                         EXCIT=(loadr),
                         OBSERVATION=_F(
                             GROUP_MA=obs_groups,
                             NOM_CHAM='DEPL',
                             NOM_CMP=('DX', 'DY', 'DZ'),
                             INST=k+dt,
                             OBSE_ETAT_INIT='NON'
                         ),
                         SCHEMA_TEMPS=_F(FORMULATION='DEPLACEMENT',
                                         SCHEMA="HHT",
                                         ALPHA=-0.3
                                         ),
                         # add damping stablize the oscilations Need to study in the future
                         INCREMENT=_F(LIST_INST=times, INST_FIN=(1+k)*dt),
                         MODELE=model,
                         )
t2 = time.time()
print(f"[TIME] DYNA_NON_LINE         : {t2 - t1:.6f} sec")


t3 = time.time()
tblp = POST_RELEVE_T(ACTION=(_F(OPERATION='EXTRACTION',      # For Extraction of values
                                INTITULE='Nodal Displacements',    # Name of the table in .resu file
                                # The result from which values will be extracted(STAT_NON_LINE)
                                RESULTAT=resn,
                                # Field to extract. DEPL = Displacements
                                NOM_CHAM=('DEPL'),
                                # TOUT_CMP='OUI',
                                # Components of DISP to extract
                                NOM_CMP=('DX', 'DY', 'DZ'), 
                                GROUP_NO='ALLNODE',               # Extract only for nodes of group DISP
                                # STAT_NON_LINE calculates for 10 INST. I want only last INST
                                INST=(1+k)*dt,
                                ),),
                     )
t4 = time.time()
print(f"[TIME] tblp    : {t4 - t3:.6f} sec")

t5 = time.time()
tblp2 = POST_RELEVE_T(ACTION=(_F(OPERATION='EXTRACTION',      # For Extraction of values
                                 INTITULE='Nodal Velocity',    # Name of the table in .resu file
                                 # The result from which values will be extracted(STAT_NON_LINE)
                                 RESULTAT=resn,
                                 # Field to extract. VITE = velocity,
                                 NOM_CHAM=('VITE'),
                                 # TOUT_CMP='OUI',
                                 # Components of DISP to extract
                                 NOM_CMP=('DX', 'DY', 'DZ'),
                                 GROUP_NO='ALLNODE',               # Extract only for nodes of group DISP
                                 # STAT_NON_LINE calculates for 10 INST. I want only last INST
                                 INST=(1+k)*dt,
                                 ),),
                      )

t6 = time.time()
print(f"[TIME] tblp2    : {t6 - t5:.6f} sec")

t7 = time.time()
posi = aster_module.get_position_aster(tblp)

#if k < itimes-1:
#    del Fnh

timeFE = dt * k

if k == 0:
    velo_nodes = np.zeros_like(posi)
else:
    velo_nodes = aster_module.get_velocity_aster(tblp2)
t8 = time.time()
print(f"[TIME] aster_module         : {t8 - t7:.6f} sec")


t9 = time.time()

# -----------------------------
# 1) wave velocity 계산 (노드 기준)
# -----------------------------
if wave_model is not None:
    wave_vel_nodes = wave_model.velocity(posi, timeFE)
    elev = wave_model.eta(posi, timeFE)
else:
    wave_vel_nodes = np.zeros_like(posi)
    elev = np.zeros(posi.shape[0], dtype=float)

# -----------------------------
# 2) current + wave
# -----------------------------
current_nodes = np.ones_like(posi, dtype=float) * np.asarray(current_input, dtype=float)
fluid_vel_nodes = current_nodes + wave_vel_nodes

# -----------------------------
# 3) net 전용: 요소 평균 유속
# -----------------------------
u_net = np.array([
    fluid_vel_nodes[[a, b, c]].mean(axis=0) if c != -1
    else fluid_vel_nodes[[a, b]].mean(axis=0)
    for (a, b, c) in netting.triangular_elements
], dtype=float)

if k == 0:
    netting.init_wake(posi, current_input)
    netting.update_wake_reduction(posi, velo_nodes)
else:
    netting.update_wake_reduction(posi, velo_nodes)

t10 = time.time()
print(f"[TIME] wake update           : {t10 - t9:.6f} sec")

t11 = time.time()

# ScreenModel: 요소 유속(u_net) 사용
net_dyn  = netting.force_on_element(posi, u_net, velo_nodes, elev)
net_buoy = netting.cal_buoy_force(posi, elev)

# MorisonModel: 노드 유속(fluid_vel_nodes) 그대로 전달
flo_dyn   = floater_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
flo_buoy  = floater_morison.cal_buoy_force(posi, elev)

btr_dyn   = bottomring_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
btr_buoy  = bottomring_morison.cal_buoy_force(posi, elev)

siro_dyn  = siderope_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
siro_buoy = siderope_morison.cal_buoy_force(posi, elev)

bri_dyn   = bridle_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
bri_buoy  = bridle_morison.cal_buoy_force(posi, elev)

buoyl_dyn  = buoyline_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
buoyl_buoy = buoyline_morison.cal_buoy_force(posi, elev)

ancla_dyn  = anchora_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
ancla_buoy = anchora_morison.cal_buoy_force(posi, elev)

anclb_dyn  = anchorb_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
anclb_buoy = anchorb_morison.cal_buoy_force(posi, elev)

moor_dyn  = mooring_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
moor_buoy = mooring_morison.cal_buoy_force(posi, elev)

bra_dyn  = braket_morison.force_on_element(posi, fluid_vel_nodes, velo_nodes, elev)
bra_buoy = braket_morison.cal_buoy_force(posi, elev)

t12 = time.time()
print(f"[TIME] force_on_element      : {t12 - t11:.6f} sec")

net_total = net_dyn + net_buoy
flo_total = flo_dyn + flo_buoy
btr_total = btr_dyn + btr_buoy
siro_total = siro_dyn + siro_buoy
bri_total = bri_dyn + bri_buoy
buoyl_total = buoyl_dyn + buoyl_buoy
ancla_total = ancla_dyn + ancla_buoy
anclb_total = anclb_dyn + anclb_buoy
moor_total = moor_dyn + moor_buoy
bra_total = bra_dyn + bra_buoy

floater_morison.hydro_total_forces = flo_total
netting.hydro_total_forces = net_total
bottomring_morison.hydro_total_forces = btr_total
siderope_morison.hydro_total_forces = siro_total
buoyline_morison.hydro_total_forces = buoyl_total
anchora_morison.hydro_total_forces = ancla_total
anchorb_morison.hydro_total_forces = anclb_total
mooring_morison.hydro_total_forces = moor_total
braket_morison.hydro_total_forces = bra_total
bridle_morison.hydro_total_forces = bri_total

outdir = os.path.join(cwd, 'pythonOutput')
os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/Fnh'), exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/posi'), exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/velo'), exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/N'), exist_ok=True)

t13 = time.time()

Fnh = (netting.distribute_force(nodenumber) 
       + floater_morison.distribute_force(nodenumber) 
       + bottomring_morison.distribute_force(nodenumber) 
       + siderope_morison.distribute_force(nodenumber) 
       + buoyline_morison.distribute_force(nodenumber)
       + anchora_morison.distribute_force(nodenumber)
       + anchorb_morison.distribute_force(nodenumber)
       + mooring_morison.distribute_force(nodenumber)
       + braket_morison.distribute_force(nodenumber)
       + bridle_morison.distribute_force(nodenumber)
       )

t14 = time.time()
print(f"[TIME] distribute_force      : {t14 - t13:.6f} sec")



t15 = time.time()

if k % save_step == 0:

    np.savetxt(os.path.join(cwd, 'pythonOutput/Fnh/Fnh_' +
                            str(round((k)*dt, 3))+'.txt'), Fnh, fmt='%.3e')
    np.savetxt(os.path.join(cwd, 'pythonOutput/posi/posi_' +
                            str(round((k)*dt, 3))+'.txt'), posi, fmt='%.3e')
    np.savetxt(os.path.join(cwd, 'pythonOutput/velo/velo_' +
                            str(round((k)*dt, 3))+'.txt'), velo_nodes, fmt='%.3e')
    

if k % save_step == 0:

    inst_now = (k+1) * dt

    grp_calc = ('NETTWIN','FLOATER','BRACKET','BOTRING', 'MOORFRM','BRIDLE','SIDEROPE','BUOYLN1','ANCHL1','ANCHL2')


    stat1 = CALC_CHAMP(
        RESULTAT=resn,
        GROUP_MA=grp_calc,
        CONTRAINTE=('EFGE_ELGA', 'EFGE_ELNO'),
    )


#     # stat_load = CALC_CHAMP(
#     #     RESULTAT=resn,
#     #     FORCE='FORC_NODA', 
#     #     INST=inst_now,
#     # )

#     # stat3 = CALC_CHAMP(
#     #     RESULTAT=resn,
#     #     FORCE='REAC_NODA',
#     # )

    reac1 = POST_RELEVE_T(
        ACTION=_F(
            GROUP_NO='ALLNODE',
            INTITULE='Element Force',
            NOM_CHAM=('EFGE_ELNO',),
            NOM_CMP=('N'), 
            OPERATION=('EXTRACTION',),
            INST=inst_now,     
            RESULTAT=stat1,
        )
    )  

    # reac3 = POST_RELEVE_T(
    #     ACTION=_F(
    #         GROUP_NO=('FLONODE',),  # floater node  FLONODE
    #         INTITULE='Reactions Forces',
    #         NOM_CHAM=('REAC_NODA',),
    #         NOM_CMP=('DX','DY','DZ'), 
    #         OPERATION=('EXTRACTION',),
    #         INST=inst_now,     
    #         RESULTAT=stat3,
    #     )
    # )         

    # node_loads = POST_RELEVE_T(
    #     ACTION=_F(
    #         OPERATION='EXTRACTION',
    #         INTITULE='Nodal Loads',
    #         RESULTAT=stat_load,
    #         NOM_CHAM='FORC_NODA', 
    #         NOM_CMP=('DX', 'DY', 'DZ'), 
    #         GROUP_NO='ALLNODE',
    #         INST=inst_now,     
    #     )
    # )                

    n = aster_module.get_N_aster(reac1)

    # # n = aster_module.get_N_aster(node_loads)

    # # IMPR_TABLE(FORMAT_R='1PE12.3',
    # #         TABLE=reac3,
    # #         UNITE=10,
    # # )

    np.savetxt(os.path.join(cwd, 'pythonOutput/N/N_' +
                            str(round((k)*dt, 3))+'.txt'), n, fmt='%.3e')

t16 = time.time()
print(f"[TIME] save file             : {t16 - t15:.6f} sec")

#print("save results.....at "+str(timeFE))

t17 = time.time()

# DETRUIRE(CONCEPT=_F(NOM=(tblp)))
# DETRUIRE(CONCEPT=_F(NOM=(tblp2)))
# if k != 0:
#     if k < itimes-1:
#         for i in range(1, nodenumber+1):
#             DETRUIRE(CONCEPT=_F(NOM=load_dynamic))

if k % 1000 == 0 and k > 0 :
    DETRUIRE(CONCEPT=_F(NOM=load_dynamic))
t18 = time.time()

print(f"[TIME] DETRUIRE file             : {t18 - t17:.6f} sec")

print(f"[TIME] TOTAL step time       : {t18 - t0:.6f} sec")
print("--------------------------------------------------")