/*******************************************************************************
* RT-Flux-PIHM is a finite volume based, reactive transport module that operates
* on top of the hydrological land surface processes described by Flux-PIHM.
* RT-Flux-PIHM tracks the transportation and reaction in a given watershed. It
* uses operator splitting technique to couple transport and reaction.
*****************************************************************************/
#include "pihm.h"

/* Begin global variable definition (MACRO) */
#define ZERO   1E-20
#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80
#define INFTYSMALL  1E-6

void InitChem(const char cdbs_filen[], const char cini_filen[],
    const ctrl_struct *ctrl, const calib_struct *cal, forc_struct *forc,
    chemtbl_struct chemtbl[], kintbl_struct kintbl[], rttbl_struct *rttbl,
    elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i, j, k;
    int             chem_ind;
    FILE           *fp;

    fp = fopen(cdbs_filen, "r");
    CheckFile(fp, cdbs_filen);

    ReadCini(cini_filen, chemtbl, rttbl->NumStc, elem);

    /*
     * Look up database to find required parameters and dependencies for
     * chemical species
     */
    Lookup(fp, cal, chemtbl, kintbl, rttbl);

    /*
     * Apply calibration
     */
    rttbl->pumps[0].Injection_rate *= cal->gwinflux;
    rttbl->pumps[0].flow_rate *= cal->gwinflux;

    for (i = 0; i < nelem; i++)
    {
        for (k = 0; k < rttbl->NumStc; k++)
        {
            elem[i].restart_input.ssa_unsat[k] *= cal->ssa;
            elem[i].restart_input.ssa_gw[k] *= cal->ssa;
        }
    }

    chem_ind = FindChem("'DOC'", chemtbl, rttbl->NumStc);
    if (chem_ind >= 0)
    {
        for (i = 0; i < nelem; i++)
        {
            elem[i].restart_input.tconc_unsat[chem_ind] *= cal->initconc;
            elem[i].restart_input.tconc_gw[chem_ind] *= cal->initconc;
        }

        rttbl->prcp_conc[chem_ind] *= cal->prcpconc;
        if (ctrl->PrpFlg == 2)
        {
            for (i = 0; i < forc->TSD_prepconc.length; i++)
            {
                forc->TSD_prepconc.data[i][chem_ind] *= cal->prcpconc;
            }
        }
    }

    /*
     * Initializing element concentrations
     */
    PIHMprintf(VL_VERBOSE, "\n Initializing concentrations... \n");

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < rttbl->NumStc; j++)
        {
            if (strcmp(chemtbl[j].ChemName, "'H+'") == 0)
            {
                elem[i].chms_unsat.t_conc[j] = elem[i].restart_input.tconc_unsat[j];
                elem[i].chms_unsat.p_actv[j] = elem[i].chms_unsat.t_conc[j];
                elem[i].chms_unsat.p_conc[j] = elem[i].chms_unsat.t_conc[j];
                elem[i].chms_unsat.ssa[j] = elem[i].restart_input.ssa_unsat[j];

                elem[i].chms_gw.t_conc[j] = elem[i].restart_input.tconc_gw[j];
                elem[i].chms_gw.p_actv[j] = elem[i].chms_gw.t_conc[j];
                elem[i].chms_gw.p_conc[j] = elem[i].chms_gw.t_conc[j];
                elem[i].chms_gw.ssa[j] = elem[i].restart_input.ssa_gw[j];
            }
            else if (chemtbl[j].itype == MINERAL)
            {
                elem[i].chms_unsat.t_conc[j] = elem[i].restart_input.tconc_unsat[j];
                /* Update the concentration of mineral using molar volume */
                elem[i].chms_unsat.t_conc[j] *= (rttbl->RelMin == 0) ?
                    /* Absolute mineral volume fraction */
                    1000.0 / chemtbl[j].MolarVolume / elem[i].soil.smcmax :
                    /* Relative mineral volume fraction */
                    (1.0 - elem[i].soil.smcmax + INFTYSMALL) * 1000.0 /
                    chemtbl[j].MolarVolume / elem[i].soil.smcmax;
                elem[i].chms_unsat.p_actv[j] = 1.0;
                elem[i].chms_unsat.p_conc[j] = elem[i].chms_unsat.t_conc[j];
                elem[i].chms_unsat.ssa[j] = elem[i].restart_input.ssa_unsat[j];

                elem[i].chms_gw.t_conc[j] = elem[i].restart_input.tconc_gw[j];
                /* Update the concentration of mineral using molar volume */
                elem[i].chms_gw.t_conc[j] *= (rttbl->RelMin == 0) ?
                    /* Absolute mineral volume fraction */
                    1000.0 / chemtbl[j].MolarVolume / elem[i].soil.smcmax :
                    /* Relative mineral volume fraction */
                    (1.0 - elem[i].soil.smcmax + INFTYSMALL) * 1000.0 /
                    chemtbl[j].MolarVolume / elem[i].soil.smcmax;
                elem[i].chms_gw.p_actv[j] = 1.0;
                elem[i].chms_gw.p_conc[j] = elem[i].chms_gw.t_conc[j];
                elem[i].chms_gw.ssa[j] = elem[i].restart_input.ssa_gw[j];
            }
            else if ((chemtbl[j].itype == CATION_ECHG) ||
                (chemtbl[j].itype == ADSORPTION))
            {
                elem[i].chms_unsat.t_conc[j] = elem[i].restart_input.tconc_unsat[j];
                elem[i].chms_unsat.p_actv[j] = elem[i].chms_unsat.t_conc[j] * 0.5;
                /* Change unit of CEC (eq g-1) into C(ion site)
                 * (eq L-1 porous space), assuming density of solid is always
                 * 2650 g L-1 */
                elem[i].chms_unsat.t_conc[j] *= (1.0 - elem[i].soil.smcmax) * 2650.0;
                elem[i].chms_unsat.p_conc[j] = elem[i].chms_unsat.t_conc[j];

                elem[i].chms_gw.t_conc[j] = elem[i].restart_input.tconc_gw[j];
                elem[i].chms_gw.p_actv[j] = elem[i].chms_gw.t_conc[j] * 0.5;
                /* Change unit of CEC (eq g-1) into C(ion site)
                 * (eq L-1 porous space), assuming density of solid is always
                 * 2650 g L-1 */
                elem[i].chms_gw.t_conc[j] *= (1.0 - elem[i].soil.smcmax) * 2650.0;
                elem[i].chms_gw.p_conc[j] = elem[i].chms_gw.t_conc[j];
            }
            else
            {
                elem[i].chms_unsat.t_conc[j] = elem[i].restart_input.tconc_unsat[j];
                elem[i].chms_unsat.p_actv[j] = elem[i].chms_unsat.t_conc[j] * 0.5;
                elem[i].chms_unsat.p_conc[j] = elem[i].chms_unsat.t_conc[j] * 0.5;
                elem[i].chms_unsat.ssa[j] = elem[i].restart_input.ssa_unsat[j];

                elem[i].chms_gw.t_conc[j] = elem[i].restart_input.tconc_gw[j];
                elem[i].chms_gw.p_actv[j] = elem[i].chms_gw.t_conc[j] * 0.5;
                elem[i].chms_gw.p_conc[j] = elem[i].chms_gw.t_conc[j] * 0.5;
                elem[i].chms_gw.ssa[j] = elem[i].restart_input.ssa_gw[j];
            }
        }

        for (j = 0; j < rttbl->NumSsc; j++)
        {
            elem[i].chms_unsat.s_conc[j] = ZERO;

            elem[i].chms_gw.s_conc[j] = ZERO;
        }
    }

    /* Speciation */
    if (!rttbl->RecFlg)
    {
        for (i = 0; i < nelem; i++)
        {
            Speciation(chemtbl, rttbl, 1, &elem[i].chms_unsat);

            Speciation(chemtbl, rttbl, 1, &elem[i].chms_gw);
        }
    }

    /* Total moles should be calculated after speciation */
    for (i = 0; i < nelem; i++)
    {
        double          vol_gw;
        double          vol_unsat;

        vol_gw = GWVol(&elem[i].topo, &elem[i].soil, &elem[i].ws);
        vol_unsat = UnsatWaterVol(&elem[i].topo, &elem[i].soil, &elem[i].ws);

        for (j = 0; j < rttbl->NumStc; j++)
        {
            if (chemtbl[j].itype == AQUEOUS)
            {
                elem[i].chms_unsat.t_mole[j] =
                    elem[i].chms_unsat.t_conc[j] * vol_unsat;

                elem[i].chms_gw.t_mole[j] =
                    elem[i].chms_gw.t_conc[j] * vol_gw;
            }
            else
            {
                elem[i].chms_unsat.t_mole[j] = 0.0;
                elem[i].chms_gw.t_mole[j] = 0.0;
            }
        }
    }

    /*
     * Initialize river concentrations
     */
    for (i = 0; i < nriver; i++)
    {
        double          vol_rivbed;
        double          vol_stream;

        vol_rivbed = RivBedVol(&river[i].topo, &river[i].matl, &river[i].ws);
        vol_stream = river[i].topo.area *
            ((river[i].ws.stage > 0.0) ? river[i].ws.stage : 0.0);

        for (k = 0; k < rttbl->NumStc; k++)
        {
            if (chemtbl[k].itype == AQUEOUS)
            {
                river[i].chms_stream.t_conc[k] =
                    0.5 * elem[river[i].leftele - 1].chms_gw.t_conc[k] +
                    0.5 * elem[river[i].rightele - 1].chms_gw.t_conc[k];
                river[i].chms_stream.p_actv[k] = river[i].chms_stream.t_conc[k];
                river[i].chms_stream.p_conc[k] = river[i].chms_stream.t_conc[k];
                river[i].chms_stream.t_mole[k] =
                    river[i].chms_stream.t_conc[k] * vol_stream;

                river[i].chms_rivbed.t_conc[k] =
                    0.5 * elem[river[i].leftele - 1].chms_gw.t_conc[k] +
                    0.5 * elem[river[i].rightele - 1].chms_gw.t_conc[k];
                river[i].chms_rivbed.p_actv[k] = river[i].chms_rivbed.t_conc[k];
                river[i].chms_rivbed.p_conc[k] = river[i].chms_rivbed.t_conc[k];
                river[i].chms_rivbed.t_mole[k] =
                    river[i].chms_rivbed.t_conc[k] * vol_rivbed;
            }
            else
            {
                river[i].chms_stream.t_conc[k] = 1.0E-20;
                river[i].chms_stream.p_conc[k] = 1.0E-20;
                river[i].chms_stream.p_actv[k] = 1.0E-20;
                river[i].chms_stream.t_mole[k] = 0.0;

                river[i].chms_rivbed.t_conc[k] = 1.0E-20;
                river[i].chms_rivbed.p_conc[k] = 1.0E-20;
                river[i].chms_rivbed.p_actv[k] = 1.0E-20;
                river[i].chms_rivbed.t_mole[k] = 0.0;
            }
        }

        for (j = 0; j < rttbl->NumSsc; j++)
        {
            river[i].chms_stream.s_conc[j] = ZERO;

            river[i].chms_rivbed.s_conc[j] = ZERO;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        for (k = 0; k < NumSpc; k++)
        {
            NV_Ith(CV_Y, UNSAT_MOLE(i, k)) = elem[i].chms_unsat.t_mole[k];
            NV_Ith(CV_Y, GW_MOLE(i, k)) = elem[i].chms_gw.t_mole[k];
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        for (k = 0; k < NumSpc; k++)
        {
            NV_Ith(CV_Y, STREAM_MOLE(i, k)) = river[i].chms_stream.t_mole[k];
            NV_Ith(CV_Y, RIVBED_MOLE(i, k)) = river[i].chms_rivbed.t_mole[k];
        }
    }

    fclose(fp);
}

void FreeChem(Chem_Data CD)
{
//    int             i;
//
//    free(CD->BTC_loc);
//    free(CD->prepconcindex);
//
//    // CD->Vcele
//    for (i = 0; i < CD->NumVol; i++)
//    {
//        free(CD->Vcele[i].log10_pconc);
//        free(CD->Vcele[i].log10_sconc);
//        free(CD->Vcele[i].p_para);
//        free(CD->Vcele[i].btcv_pconc);
//    }
//    free(CD->Vcele);
//
//    free(CD->Flux);
//
//    if (CD->NumPUMP > 0)
//    {
//        free(CD->pumps);
//    }
//
//    // CD->TSD_prepconc
//    for (i = 0; i < CD->TSD_prepconc[0].length; i++)
//    {
//        free(CD->TSD_prepconc[0].data[i]);
//    }
//    free(CD->TSD_prepconc[0].data);
//    free(CD->TSD_prepconc[0].ftime);
//    free(CD->TSD_prepconc[0].value);
//    free(CD->TSD_prepconc);
//
//    free(CD->Precipitation.chms.t_conc);
//    free(CD->Precipitation.chms.p_conc);
//    free(CD->Precipitation.p_para);
//
}

double GWVol(const topo_struct *topo, const soil_struct *soil,
    const wstate_struct *ws)
{
    double          vol;

    if (ws->gw < 0.0)
    {
        vol = 0.0;
    }
    else if (ws->gw > soil->depth)
    {
        vol = soil->depth * soil->smcmax + (ws->gw - soil->depth) * soil->porosity;
    }
    else
    {
        vol = ws->gw * soil->smcmax;
    }

    return vol * topo->area;
}

double UnsatWaterVol(const topo_struct *topo, const soil_struct *soil,
    const wstate_struct *ws)
{
    double          deficit;
    double          vol;

    deficit = soil->depth - ws->gw;
    deficit = (deficit < soil->depth) ? deficit : soil->depth;
    deficit = (deficit > 0.0) ? deficit : 0.0;

    vol = deficit * soil->smcmin +
        ((ws->unsat < 0.0) ? 0.0 : ws->unsat * soil->porosity);

    return vol * topo->area;
}

double UnsatSatRatio(double depth, double unsat, double gw)
{
    return ((unsat < 0.0) ? 0.0 : ((gw > depth) ? 1.0 : unsat / (depth - gw)));
}

double RivBedVol(const river_topo_struct *topo, const matl_struct *matl,
    const river_wstate_struct *ws)
{
    double          vol;

    if (ws->gw < 0.0)
    {
        vol = 0.0;
    }
    else if (ws->gw > matl->bedthick)
    {
        vol = matl->bedthick * (matl->porosity + matl->smcmin) +
            (ws->gw - matl->bedthick) * matl->porosity;
    }
    else
    {
        vol = ws->gw * (matl->porosity + matl->smcmin);
    }

    return vol * topo->area;
}
