import time
t0 = time.time()

save_interval = 0.01
save_step = int(save_interval / dt)

loadr = []
loadr.append(_F(CHARGE=gF), )
loadr.append(_F(CHARGE=fixed), )
loadr.append(_F(CHARGE=sinkF1), )


t1 = time.time()
if k == 0:
    resn = DYNA_NON_LINE(
        CARA_ELEM=elemprop,
        CHAM_MATER=fieldmat,
        # ARCHIVAGE=_F(PAS_ARCH=1,),
        COMPORTEMENT=(
            _F(DEFORMATION='GROT_GDEP',
               GROUP_MA=('NETTWIN',),
               RELATION='CABLE'),
            _F(DEFORMATION='PETIT',
               GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING'),
               RELATION='ELAS'),
        ),
        CONVERGENCE=_F(ITER_GLOB_MAXI=100,
                       RESI_GLOB_RELA=1e-3),
        EXCIT=(loadr),
        OBSERVATION=_F(GROUP_MA=('NETTWIN','FLOATER', 'BRACKET', 'BOTRING'),
                        NOM_CHAM='DEPL',
                        NOM_CMP=('DX', 'DY', 'DZ'),
                        INST=k+dt,
                        OBSE_ETAT_INIT='NON'),
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
    loadr.append(_F(CHARGE=load_dynamic))

    resn = DYNA_NON_LINE(CARA_ELEM=elemprop,
                         CHAM_MATER=fieldmat,
                         ETAT_INIT=_F(EVOL_NOLI=resn),
                        #  ARCHIVAGE=_F( PAS_ARCH=10,),
                         COMPORTEMENT=(_F(DEFORMATION='GROT_GDEP',
                                          GROUP_MA=('NETTWIN',),
                                          RELATION='CABLE'),
                                       _F(DEFORMATION='PETIT',
                                          GROUP_MA=('FLOATER', 'BRACKET', 'BOTRING'),
                                          RELATION='ELAS'),
                                       ),
                         CONVERGENCE=_F(ITER_GLOB_MAXI=100,
                                        RESI_GLOB_RELA=1e-3),
                         EXCIT=(loadr),
                        OBSERVATION=_F(GROUP_MA=('NETTWIN','FLOATER', 'BRACKET', 'BOTRING'),
                                        NOM_CHAM='DEPL',
                                         NOM_CMP=('DX', 'DY', 'DZ'),
                                         INST=k+dt,
                                         OBSE_ETAT_INIT='NON'),
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

Nnode = posi.shape[0]
elev = np.zeros(Nnode)

if k == 0:
    velo_nodes = np.zeros_like(posi)
else:
    velo_nodes = aster_module.get_velocity_aster(tblp2)
t8 = time.time()
print(f"[TIME] aster_module         : {t8 - t7:.6f} sec")


t9 = time.time()
u = np.ones((NO_net, 3)) * current_input

if k == 0:
    netting.init_wake(posi, current_input) 
    netting.update_wake_reduction(posi, velo_nodes)
else:
    netting.update_wake_reduction(posi, velo_nodes)
t10 = time.time()
print(f"[TIME] wake update           : {t10 - t9:.6f} sec")

t11 = time.time()

f_dyn  = netting.force_on_element(posi, u, velo_nodes, elev)
f_buoy = netting.cal_buoy_force(posi, elev)

t12 = time.time()
print(f"[TIME] force_on_element      : {t12 - t11:.6f} sec")

f_total = f_dyn + f_buoy
netting.hydro_total_forces = f_total

outdir = os.path.join(cwd, 'pythonOutput')
os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/Fnh'), exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/posi'), exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/velo'), exist_ok=True)
os.makedirs(os.path.join(cwd, 'pythonOutput/N'), exist_ok=True)

t13 = time.time()

Fnh = netting.distribute_force(nodenumber)

t14 = time.time()
print(f"[TIME] distribute_force      : {t14 - t13:.6f} sec")



t15 = time.time()

if k % save_step == 0:

    np.savetxt(
        os.path.join(outdir, f'netting_{round(k*dt,3)}.txt'),
        f_total,
        fmt='%.3e'
    )
    np.savetxt(os.path.join(cwd, 'pythonOutput/Fnh/Fnh_' +
                            str(round((k)*dt, 3))+'.txt'), Fnh, fmt='%.3e')
    np.savetxt(os.path.join(cwd, 'pythonOutput/posi/posi_' +
                            str(round((k)*dt, 3))+'.txt'), posi, fmt='%.3e')
    np.savetxt(os.path.join(cwd, 'pythonOutput/velo/velo_' +
                            str(round((k)*dt, 3))+'.txt'), velo_nodes, fmt='%.3e')
    

if k % save_step == 0:

    inst_now = (k+1) * dt

    stat1 = CALC_CHAMP(
        RESULTAT=resn,
        GROUP_MA=('NETTWIN','FLOATER', 'BRACKET', 'BOTRING',),
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