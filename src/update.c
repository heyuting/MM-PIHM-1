
/*****************************************************************************
 * File		: update.c
 * Function	: Update variables during simulation
 ****************************************************************************/

#include "pihm.h"

void Summary (pihm_struct pihm, N_Vector CV_Y, double stepsize)
{
    double         *y;
    double          wtd0, wtd1, elemsatn0, elemsatn1, realunsat0, realunsat1,
        realgw0, realgw1, recharge, runoff;
    double          totalw0, totalw1;
    int             i, j;
    elem_struct    *elem;

    y = NV_DATA_S (CV_Y);


    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        elem->unsat = y[i + pihm->numele];
        elem->gw = y[i + 2 * pihm->numele];

        /* calculate infiltration based on mass conservation */
        wtd0 = elem->soil.depth - ((elem->gw0 > 0.0) ? elem->gw0 : 0.0);
        wtd0 = (wtd0 < 0.0) ? 0.0 : wtd0;
        if (wtd0 <= 0.0)
        {
            elemsatn0 = 1.0;
        }
        else if (elem->unsat0 < 0.0)
        {
            elemsatn0 = 0.0;
        }
        else
        {
            elemsatn0 = elem->unsat0 / wtd0;
        }
        elemsatn0 = (elemsatn0 > 1.0) ? 1.0 : elemsatn0;
        elemsatn0 = (elemsatn0 < 0.0) ? 0.0 : elemsatn0;
        realunsat0 = elemsatn0 * wtd0;

        realgw0 = elem->gw0;
        realgw0 = (realgw0 > elem->soil.depth) ? elem->soil.depth : realgw0;
        realgw0 = (realgw0 < 0.0) ? 0.0 : realgw0;

        totalw0 = elem->gw0 + elem->unsat0;
        totalw0 = (totalw0 > elem->soil.depth) ? elem->soil.depth : totalw0;
        totalw0 = (totalw0 < 0.0) ? 0.0 : totalw0;

        wtd1 = elem->soil.depth - ((elem->gw > 0.0) ? elem->gw : 0.0);
        wtd1 = (wtd1 < 0.0) ? 0.0 : wtd1;
        if (wtd1 <= 0.0)
        {
            elemsatn1 = 1.0;
        }
        else if (elem->unsat < 0.0)
        {
            elemsatn1 = 0.0;
        }
        else
        {
            elemsatn1 = elem->unsat / wtd1;
        }
        elemsatn1 = (elemsatn1 > 1.0) ? 1.0 : elemsatn1;
        elemsatn1 = (elemsatn1 > 0.0) ? 0.0 : elemsatn1;
        realunsat1 = elemsatn1 * wtd1;

        realgw1 = (elem->gw > elem->soil.depth) ? elem->soil.depth : elem->gw;

        totalw1 = elem->gw + elem->unsat;
        totalw1 = (totalw1 > elem->soil.depth) ? elem->soil.depth : totalw1;
        totalw1 = (totalw1 < 0.0) ? 0.0 : totalw1;
        /* subsurface runoff rate */
        elem->runoff = 0.0;
        for (j = 0; j < 3; j++)
        {
            elem->runoff += elem->fluxsub[j] / elem->topo.area;
        }
#ifdef _NOAH_
        recharge =
            (realgw1 - realgw0) * elem->soil.porosity / stepsize +
            elem->runoff + elem->et_from_sat * elem->et[1];
        //elem->infil = (realunsat1 - realunsat0) * elem->soil.porosity / stepsize + recharge + (1.0 - elem->et_from_sat) * elem->et[1] + elem->et[2];
        elem->infil =
            (totalw1 - totalw0) * elem->soil.porosity / stepsize +
            elem->runoff + elem->et[1] + elem->et[2];
#else
        recharge =
            (realgw1 - realgw0) * elem->soil.porosity / stepsize +
            elem->runoff + ((elem->gw0 >
                aquiferdepth - elem->lc.rzd) ? elem->et[1] : 0.0);
        elem->infil =
            (realunsat1 - realunsat0) * elem->soil.porosity / stepsize +
            recharge + (elem->surf0 <
            EPS / 100.0 ? elem->et[2] : 0.0) + ((elem->gw0 <=
                aquiferdepth - elem->lc.rzd) ? elem->et[1] : 0.0);
#endif
        if (elem->infil < 0.0)
        {
            elem->runoff -= elem->infil;
            elem->infil = 0.0;
        }

        for (j = 0; j < 3; j++)
        {
            elem->fluxtotal[j] = elem->fluxsurf[j] + elem->fluxsub[j];
        }
    }
    for (i = 0; i < pihm->numele; i++)
    {
        pihm->elem[i].surf0 = y[i];
        pihm->elem[i].unsat0 = y[i + pihm->numele];
        pihm->elem[i].gw0 = y[i + 2 * pihm->numele];
    }
    for (i = 0; i < pihm->numriv; i++)
    {
        pihm->riv[i].stage0 = y[i + 3 * pihm->numele];
        pihm->riv[i].gw0 = y[i + 3 * pihm->numele + pihm->numriv];
    }
}
