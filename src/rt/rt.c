/*******************************************************************************
* RT-Flux-PIHM is a finite volume based, reactive transport module that operates
* on top of the hydrological land surface processes described by Flux-PIHM.
* RT-Flux-PIHM tracks the transportation and reaction in a given watershed. It
* uses operator splitting technique to couple transport and reaction.
*****************************************************************************/
#include "pihm.h"

void InitChem(const char cdbs_filen[], const calib_struct *cal,
    forc_struct *forc, chemtbl_struct chemtbl[], kintbl_struct kintbl[],
    rttbl_struct *rttbl)
{
    int             i;
    int             chem_ind;
    FILE           *fp;

    fp = fopen(cdbs_filen, "r");
    CheckFile(fp, cdbs_filen);

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

    chem_ind = FindChem("'DOC'", chemtbl, rttbl->NumStc);
    if (chem_ind >= 0)
    {
        rttbl->prcp_conc[chem_ind] *= cal->prcpconc;

        if (forc->PrpFlg == 2)
        {
            for (i = 0; i < forc->TSD_prepconc.length; i++)
            {
                forc->TSD_prepconc.data[i][chem_ind] *= cal->prcpconc;
            }
        }
    }

    fclose(fp);
}

void InitRTVar(chemtbl_struct chemtbl[], rttbl_struct *rttbl,
    elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i;

    /*
     * Initializing element concentrations
     */
    PIHMprintf(VL_VERBOSE, "\n Initializing concentrations... \n");

    for (i = 0; i < nelem; i++)
    {
        int             j;

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
                    (1.0 - elem[i].soil.smcmax) * 1000.0 /
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
                    (1.0 - elem[i].soil.smcmax) * 1000.0 /
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
            elem[i].chms_unsat.s_conc[j] = ZERO_CONC;
            elem[i].chms_gw.s_conc[j] = ZERO_CONC;
        }

        /* Speciation */
        if (rttbl->RecFlg == KIN_REACTION)
        {
            _Speciation(chemtbl, rttbl, 1, &elem[i].chms_unsat);

            _Speciation(chemtbl, rttbl, 1, &elem[i].chms_gw);
        }
    }

    /* Total moles should be calculated after speciation */
    for (i = 0; i < nelem; i++)
    {
        double          vol_gw;
        double          vol_unsat;
        int             k;

        vol_gw = MAX(GWStrg(&elem[i].soil, &elem[i].ws), DEPTHR) * elem[i].topo.area;
        vol_unsat = MAX(UnsatWaterStrg(&elem[i].soil, &elem[i].ws), DEPTHR) * elem[i].topo.area;

        for (k = 0; k < rttbl->NumStc; k++)
        {
            if (chemtbl[k].itype == AQUEOUS)
            {
                elem[i].chms_unsat.t_mole[k] =
                    elem[i].chms_unsat.t_conc[k] * vol_unsat;

                elem[i].chms_gw.t_mole[k] =
                    elem[i].chms_gw.t_conc[k] * vol_gw;
            }
            else
            {
                elem[i].chms_unsat.t_mole[k] = 0.0;
                elem[i].chms_gw.t_mole[k] = 0.0;
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
        int             k;

        vol_rivbed = MAX(RivBedStrg(&river[i].matl, &river[i].ws), DEPTHR) * river[i].topo.area;
        vol_stream = river[i].topo.area * MAX(river[i].ws.stage, DEPTHR);

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
                river[i].chms_stream.t_conc[k] = ZERO_CONC;
                river[i].chms_stream.p_conc[k] = ZERO_CONC;
                river[i].chms_stream.p_actv[k] = ZERO_CONC;
                river[i].chms_stream.t_mole[k] = 0.0;

                river[i].chms_rivbed.t_conc[k] = ZERO_CONC;
                river[i].chms_rivbed.p_conc[k] = ZERO_CONC;
                river[i].chms_rivbed.p_actv[k] = ZERO_CONC;
                river[i].chms_rivbed.t_mole[k] = 0.0;
            }
        }

        for (k = 0; k < rttbl->NumSsc; k++)
        {
            river[i].chms_stream.s_conc[k] = ZERO_CONC;
            river[i].chms_rivbed.s_conc[k] = ZERO_CONC;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            NV_Ith(CV_Y, UNSAT_MOLE(i, k)) = elem[i].chms_unsat.t_mole[k];
            NV_Ith(CV_Y, GW_MOLE(i, k)) = elem[i].chms_gw.t_mole[k];

            elem[i].chmf.react_unsat[k] = 0.0;
            elem[i].chmf.react_gw[k] = 0.0;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            NV_Ith(CV_Y, STREAM_MOLE(i, k)) = river[i].chms_stream.t_mole[k];
            NV_Ith(CV_Y, RIVBED_MOLE(i, k)) = river[i].chms_rivbed.t_mole[k];
        }
    }
}

double GWStrg(const soil_struct *soil, const wstate_struct *ws)
{
    double          strg;

    if (ws->gw < 0.0)
    {
        strg = 0.0;
    }
    else if (ws->gw > soil->depth)
    {
        strg = soil->depth * soil->smcmax +
            (ws->gw - soil->depth) * soil->porosity;
    }
    else
    {
        strg = ws->gw * soil->smcmax;
    }

    return strg;
}

double UnsatWaterStrg(const soil_struct *soil, const wstate_struct *ws)
{
    double          deficit;

    deficit = soil->depth - ws->gw;
    deficit = MIN(deficit, soil->depth);
    deficit = MAX(deficit, 0.0);

    return deficit * soil->smcmin + MAX(ws->unsat, 0.0) * soil->porosity;
}

double UnsatSatRatio(double depth, double unsat, double gw)
{
    return ((unsat < 0.0) ? 0.0 : ((gw > depth) ? 1.0 : unsat / (depth - gw)));
}

double RivBedStrg(const matl_struct *matl, const river_wstate_struct *ws)
{
    double          strg;

    if (ws->gw < 0.0)
    {
        strg = 0.0;
    }
    else if (ws->gw > matl->bedthick)
    {
        strg = matl->bedthick * (matl->porosity + matl->smcmin) +
            (ws->gw - matl->bedthick) * matl->porosity;
    }
    else
    {
        strg = ws->gw * (matl->porosity + matl->smcmin);
    }

    return strg;
}
