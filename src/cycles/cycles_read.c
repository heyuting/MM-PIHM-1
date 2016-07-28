#include "pihm.h"

void ReadCyclesCtrl (char *filename, agtbl_struct *agtbl, ctrl_struct *ctrl,
    int numele)
{
    FILE           *simctrl_file;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             i;
    int             match;
    int             index;

    /* Open simulation control file */
    simctrl_file = fopen (filename, "r");
    CheckFile (simctrl_file, filename);

    agtbl->op = (int *)malloc (numele * sizeof (int));
    agtbl->rotsz = (int *)malloc (numele * sizeof (int));
    agtbl->auto_N = (int *)malloc (numele * sizeof (int));
    agtbl->auto_P = (int *)malloc (numele * sizeof (int));
    agtbl->auto_S = (int *)malloc (numele * sizeof (int));

    /* Read simulation control file */
    FindLine (simctrl_file, "BOF");

    NextLine (simctrl_file, cmdstr);
    for (i = 0; i < numele; i++)
    {
        NextLine (simctrl_file, cmdstr);
        match =
            sscanf (cmdstr, "%d %d %d %d %d %d", &index, &agtbl->op[i],
            &agtbl->rotsz[i], &agtbl->auto_N[i], &agtbl->auto_P[i],
            &agtbl->auto_S[i]);
        if (match != 6)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".cycles file format error!\n");
            exit (1);
        }
    }

    FindLine (simctrl_file, "OPERATION_FILE");

    i = 0;
    while (1)
    {
        NextLine (simctrl_file, cmdstr);
        sscanf (cmdstr, "%s", optstr);

        if (strcasecmp (cmdstr, "EOF") == 0 ||
            strcasecmp (optstr, "PRINT_CTRL") == 0)
        {
            break;
        }

        match = sscanf (cmdstr, "%d %s", &index, agtbl->opfilen[i]);
        if (match != 2 || i != index - 1)
        {
            fprintf (stderr, "Error reading %s.\n", filename),
                fprintf (stderr, "Please check file format.\n");
            PIHMExit (EXIT_FAILURE);
        }
        i++;
    }

    agtbl->nopfile = i;

    /* Output control */
    FindLine (simctrl_file, "PRINT_CTRL");

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "BIOMASS", &ctrl->prtvrbl[BIOMASS_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "LAI", &ctrl->prtvrbl[LAI_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "RADN_INTCP", &ctrl->prtvrbl[RADNINTCP_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "WATER_STRESS", &ctrl->prtvrbl[WATER_STS_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "N_STRESS", &ctrl->prtvrbl[N_STS_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "CROP_TR", &ctrl->prtvrbl[CROP_TR_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "CROP_POT_TR", &ctrl->prtvrbl[CROP_POTTR_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "RES_EVAP", &ctrl->prtvrbl[RES_EVAP_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NO3_PROF", &ctrl->prtvrbl[NO3_PROF_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NO3_RIVER", &ctrl->prtvrbl[NO3_RIVER_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NH4_PROF", &ctrl->prtvrbl[NH4_PROF_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NH4_RIVER", &ctrl->prtvrbl[NH4_RIVER_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NO3_DENITRIF", &ctrl->prtvrbl[NO3_DENIT_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NO3_LEACH", &ctrl->prtvrbl[NO3_LEACH_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NH4_LEACH", &ctrl->prtvrbl[NH4_LEACH_CTRL],
            'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NO3_LEACH_RIVER",
        &ctrl->prtvrbl[NO3_LEACH_RIVER_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (simctrl_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NH4_LEACH_RIVER",
        &ctrl->prtvrbl[NH4_LEACH_RIVER_CTRL], 'i'))
    {
        fprintf (stderr, "Please check file format.\n");
        PIHMExit (EXIT_FAILURE);
    }

    fclose (simctrl_file);
}

void ReadSoilInit (char *filename, soiltbl_struct *soiltbl)
{
    FILE           *soil_file;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             layer;
    int             i, j;

    /* 
     * Open soil initialization file
     */
    soil_file = fopen (filename, "r");
    CheckFile (soil_file, filename);

    soiltbl->totalLayers = (int *)malloc (soiltbl->number * sizeof (int));
    soiltbl->clay_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->sand_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->iom_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->bd_lyr = (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->NO3_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->NH4_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));

    /* Read soil file */
    FindLine (soil_file, "BOF");

    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine (soil_file, cmdstr);
        if (!ReadKeyword (cmdstr, "SOIL_TYPE", &index, 'i'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        if (i != index - 1)
        {
            printf ("Cannot read information of the %dth soil type!\n",
                i + 1);
            printf (".soilinit file format error!\n");
            exit (1);
        }

        NextLine (soil_file, cmdstr);
        if (!ReadKeyword (cmdstr, "TOTAL_LAYERS", &soiltbl->totalLayers[i],
                'i'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        soiltbl->clay_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->sand_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->iom_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->bd_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->NO3_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->NH4_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));

        /* Skip header */
        NextLine (soil_file, cmdstr);

        for (j = 0; j < soiltbl->totalLayers[i]; j++)
        {
            NextLine (soil_file, cmdstr);
            match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf", &layer,
                &soiltbl->clay_lyr[i][j], &soiltbl->sand_lyr[i][j],
                &soiltbl->iom_lyr[i][j], &soiltbl->bd_lyr[i][j],
                &soiltbl->NO3_lyr[i][j], &soiltbl->NH4_lyr[i][j]);

            if (match != 7 || j != layer - 1)
            {
                printf
                    ("Cannot read information of the %dth layer of the %dth"
                    "soil type!\n", j + 1, i + 1);
                printf (".soilinit file format error!\n");
                exit (1);
            }
        }
    }

    fclose (soil_file);
}

void ReadCrop (char *filename, croptbl_struct *croptbl)
{
    FILE           *crop_file;
    char            cmdstr[MAXSTRING];
    char            temp[MAXSTRING];
    int             j;

    crop_file = fopen (filename, "r");
    CheckFile (crop_file, filename);

    /* Read crop description file */
    /* First count how many crop types are there in the description file */
    croptbl->number = CountOccurance (crop_file, "NAME");

    croptbl->cropName = (char **)malloc (croptbl->number * sizeof (char *));
    croptbl->userFloweringTT =
        (double *)malloc (croptbl->number * sizeof (int));
    croptbl->userMaturityTT =
        (double *)malloc (croptbl->number * sizeof (int));
    croptbl->userMaximumSoilCoverage =
        (double *)malloc (croptbl->number * sizeof (double));

    croptbl->userMaximumRootingDepth =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userExpectedYieldAvg =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userExpectedYieldMax =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userExpectedYieldMin =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userPercentMoistureInYield =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userFractionResidueStanding =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userFractionResidueRemoved =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingBiomassThresholdUpper =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingBiomassThresholdLower =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingTiming =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingDestiny =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userTranspirationMinTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTranspirationThresholdTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userColdDamageMinTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userColdDamageThresholdTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTemperatureBase =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTemperatureOptimum =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTemperatureMaximum =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userShootPartitionInitial =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userShootPartitionFinal =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userRadiationUseEfficiency =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTranspirationUseEfficiency =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIx = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIo = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIk = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userEmergenceTT =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userNMaxConcentration =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userNDilutionSlope =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userKc = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userAnnual = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userLegume = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userC3orC4 = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userExtinctionCoefficient =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userPlantingDensity =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingStart =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userClippingEnd = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->LWP_StressOnset =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->LWP_WiltingPoint =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->transpirationMax =
        (double *)malloc (croptbl->number * sizeof (double));

    /* Rewind to the beginning of file */
    FindLine (crop_file, "BOF");

    /* Read crop properties */
    for (j = 0; j < croptbl->number; j++)
    {
        croptbl->cropName[j] = (char *)malloc (MAXSTRING * sizeof (char));
        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "NAME", croptbl->cropName[j], 's'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "FLOWERING_TT",
                &croptbl->userFloweringTT[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MATURITY_TT", &croptbl->userMaturityTT[j],
                'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MAXIMUM_SOIL_COVERAGE",
                &croptbl->userMaximumSoilCoverage[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MAXIMUM_ROOTING_DEPTH",
                &croptbl->userMaximumRootingDepth[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "AVERAGE_EXPECTED_YIELD",
                &croptbl->userExpectedYieldAvg[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MAXIMUM_EXPECTED_YIELD",
                &croptbl->userExpectedYieldMax[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MINIMUM_EXPECTED_YIELD",
                &croptbl->userExpectedYieldMin[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "COMMERCIAL_YIELD_MOISTURE",
                &croptbl->userPercentMoistureInYield[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "STANDING_RESIDUE_AT_HARVEST",
                &croptbl->userFractionResidueStanding[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "RESIDUE_REMOVED",
                &croptbl->userFractionResidueRemoved[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_UPPER",
                &croptbl->userClippingBiomassThresholdUpper[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_LOWER",
                &croptbl->userClippingBiomassThresholdLower[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "HARVEST_TIMING",
                &croptbl->userClippingTiming[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "CLIPPING_BIOMASS_DESTINY", temp, 's'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }
        if (strcasecmp ("REMOVE", temp) == 0)
        {
            croptbl->userClippingDestiny[j] = REMOVE_CLIPPING;
        }
        else if (strcasecmp ("RETURN", temp) == 0)
        {
            croptbl->userClippingDestiny[j] = RETURN_CLIPPING;
        }
        else if (strcasecmp ("GRAZING", temp) == 0)
        {
            croptbl->userClippingDestiny[j] = GRAZING_CLIPPING;
        }
        else
        {
            printf ("Option %s not recoganized!\n", temp);
            exit (1);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MIN_TEMPERATURE_FOR_TRANSPIRATION",
                &croptbl->userTranspirationMinTemperature[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "THRESHOLD_TEMPERATURE_FOR_TRANPIRATION",
                &croptbl->userTranspirationThresholdTemperature[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MIN_TEMPERATURE_FOR_COLD_DAMAGE",
                &croptbl->userColdDamageMinTemperature[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "THRESHOLD_TEMPERATURE_FOR_COLD_DAMAGE",
                &croptbl->userColdDamageThresholdTemperature[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "BASE_TEMPERATURE_FOR_DEVELOPMENT",
                &croptbl->userTemperatureBase[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "OPTIMUM_TEMPERATURE_FOR_DEVELOPEMENT",
                &croptbl->userTemperatureOptimum[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MAX_TEMPERATURE_FOR_DEVELOPMENT",
                &croptbl->userTemperatureMaximum[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "INITIAL_PARTITIONING_TO_SHOOT",
                &croptbl->userShootPartitionInitial[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "FINAL_PARTITIONING_TO_SHOOT",
                &croptbl->userShootPartitionFinal[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "RADIATION_USE_EFFICIENCY",
                &croptbl->userRadiationUseEfficiency[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "TRANSPIRATION_USE_EFFICIENCY",
                &croptbl->userTranspirationUseEfficiency[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MAXIMUM_HARVEST_INDEX",
                &croptbl->userHIx[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "MINIMUM_HARVEST_INDEX",
                &croptbl->userHIo[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "HARVEST_INDEX", &croptbl->userHIk[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "THERMAL_TIME_TO_EMERGENCE",
                &croptbl->userEmergenceTT[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "N_MAX_CONCENTRATION",
                &croptbl->userNMaxConcentration[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "N_DILUTION_SLOPE",
                &croptbl->userNDilutionSlope[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "KC", &croptbl->userKc[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "ANNUAL", &croptbl->userAnnual[j], 'i'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "LEGUME", &croptbl->userLegume[j], 'i'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "C3", &croptbl->userC3orC4[j], 'i'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "LWP_STRESS_ONSET",
                &croptbl->LWP_StressOnset[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "LWP_WILTING_POINT",
                &croptbl->LWP_WiltingPoint[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        NextLine (crop_file, cmdstr);
        if (!ReadKeyword (cmdstr, "TRANSPIRATION_MAX",
                &croptbl->transpirationMax[j], 'd'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }
    }

    fclose (crop_file);
}

void ReadOperation (const agtbl_struct *agtbl, mgmttbl_struct *mgmttbl,
    const croptbl_struct *croptbl)
{
    FILE           *op_file;
    char            cmdstr[MAXSTRING];
    char            filename[MAXSTRING];
    int             ntill;
    int             nplnt;
    int             nirrg;
    int             nfert;
    int             nautoirrg;
    int             i, j, k;
    cropmgmt_struct *cropmgmt;
    op_struct      *q;

    mgmttbl->number = agtbl->nopfile;

    mgmttbl->cropmgmt =
        (cropmgmt_struct *)malloc (mgmttbl->number *
        sizeof (cropmgmt_struct));

    for (i = 0; i < mgmttbl->number; i++)
    {
        cropmgmt = &mgmttbl->cropmgmt[i];

        cropmgmt->usingAutoIrr = 0;

        sprintf (filename, "input/%s/%s", project, agtbl->opfilen[i]);
        op_file = fopen (filename, "r");
        CheckFile (op_file, filename);

        FindLine (op_file, "BOF");
        nplnt = CountOccurance (op_file, "PLANTING");

        FindLine (op_file, "BOF");
        ntill = CountOccurance (op_file, "TILLAGE");

        FindLine (op_file, "BOF");
        nirrg = CountOccurance (op_file, "FIXED_IRRIGATION");

        FindLine (op_file, "BOF");
        nfert = CountOccurance (op_file, "FIXED_FERTILIZATION");

        FindLine (op_file, "BOF");
        nautoirrg = CountOccurance (op_file, "AUTO_IRRIGATION");

        cropmgmt->totalCropsPerRotation = nplnt;
        if (nplnt > 0)
        {
            cropmgmt->plantingOrder =
                (op_struct *)malloc (nplnt * sizeof (op_struct));

            /* Rewind to the beginning of file and read all planting operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nplnt; j++)
            {
                q = &(cropmgmt->plantingOrder[j]);

                FindLine (op_file, "PLANTING");

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "DOY", &q->opDay, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "CROP", q->cropName, 's'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "USE_AUTO_IRR",
                        &q->usesAutoIrrigation, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
                if (q->usesAutoIrrigation == 0)
                {
                    q->usesAutoIrrigation = -1;
                }
                else
                {
                    cropmgmt->usingAutoIrr = 1;
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "USE_AUTO_FERT",
                        &q->usesAutoFertilization, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
                if (q->usesAutoFertilization == 0)
                {
                    q->usesAutoFertilization = -1;
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "FRACTION", &q->plantingDensity,
                        'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "CLIPPING_START", &q->clippingStart,
                        'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
                if (q->clippingStart > 366 || q->clippingStart < 1)
                {
                    printf
                        ("ERROR: Please specify valid DOY for clipping start date!\n");
                    exit (1);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "CLIPPING_END", &q->clippingEnd,
                        'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
                if (q->clippingEnd > 366 || q->clippingEnd < 1)
                {
                    printf
                        ("ERROR: Please specify valid DOY for clipping end date!\n");
                    exit (1);
                }

                q->status = 0;

                /* Link planting order and crop description */
                for (k = 0; k < croptbl->number; k++)
                {
                    if (strcmp (q->cropName, croptbl->cropName[k]) == 0)
                    {
                        q->plantID = k;
                        break;
                    }
                }
                if (k >= croptbl->number)
                {
                    printf
                        ("ERROR: Cannot find the plant description of %s, please check your input file\n",
                        q->cropName);
                    exit (1);
                }
            }
        }

        cropmgmt->numTillage = ntill;
        if (ntill > 0)
        {
            cropmgmt->Tillage =
                (op_struct *)malloc (ntill * sizeof (op_struct));

            /* Rewind to the beginning of file and read all tillage operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < ntill; j++)
            {
                q = &(cropmgmt->Tillage[j]);

                FindLine (op_file, "TILLAGE");

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "DOY", &q->opDay, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "TOOL", q->opToolName, 's'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "DEPTH", &q->opDepth, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "SOIL_DISTURB_RATIO", &q->opSDR,
                        'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "MIXING_EFFICIENCY",
                        &q->opMixingEfficiency, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "CROP_NAME", q->cropNameT, 's'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                /* Check if the specified crop exists */
                if (strcasecmp (q->cropNameT, "N/A") != 0 &&
                    strcasecmp (q->cropNameT, "All") != 0 &&
                    !CropExist (q->cropNameT, croptbl))
                {
                    printf ("ERROR: Crop name %s not recognized!\n",
                        q->cropNameT);
                    exit (1);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "FRAC_THERMAL_TIME",
                        &q->fractionThermalTime, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "KILL_EFFICIENCY",
                        &q->killEfficiency, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "GRAIN_HARVEST", &q->grainHarvest,
                        'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "FORAGE_HARVEST",
                        &q->forageHarvest, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                q->status = 0;
            }
        }

        cropmgmt->numFertilization = nfert;
        if (nfert > 0)
        {
            cropmgmt->FixedFertilization =
                (op_struct *)malloc (nfert * sizeof (op_struct));

            /* Rewind to the beginning of file and read all fertilization
             * operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nfert; j++)
            {
                q = &(cropmgmt->FixedFertilization[j]);

                FindLine (op_file, "FIXED_FERTILIZATION");

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "DOY", &q->opDay, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "SOURCE", q->opSource, 's'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "MASS", &q->opMass, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "FORM", q->opForm, 's'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "METHOD", q->opMethod, 's'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "LAYER", &q->opLayer, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "C_ORGANIC", &q->opC_Organic, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "C_CHARCOAL", &q->opC_Charcoal,
                        'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "N_ORGANIC", &q->opN_Organic, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "N_CHARCOAL", &q->opN_Charcoal,
                        'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "N_NH4", &q->opN_NH4, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "N_NO3", &q->opN_NO3, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "P_ORGANIC", &q->opP_Organic, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "P_CHARCOAL", &q->opP_Charcoal,
                        'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "P_INORGANIC", &q->opP_Inorganic,
                        'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "K", &q->opK, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "S", &q->opS, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                q->status = 0;

                if (q->opC_Organic + q->opC_Charcoal + q->opN_Organic +
                    q->opN_Charcoal + q->opN_NH4 + q->opN_NO3 +
                    q->opP_Organic + q->opP_Charcoal + q->opP_Inorganic +
                    q->opK + q->opS <= 1.0)
                {
                    q->opMass /= 1000.0;
                }
                else
                {
                    printf
                        ("ERROR: Added fertilization fractions must be <= 1\n");
                    exit (1);
                }
            }
        }

        cropmgmt->numIrrigation = nirrg;
        if (nirrg > 0)
        {
            cropmgmt->FixedIrrigation =
                (op_struct *)malloc (nirrg * sizeof (op_struct));

            /* Rewind to the beginning of file and read all irrigation
             * operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nirrg; j++)
            {
                q = &(cropmgmt->FixedIrrigation[j]);

                FindLine (op_file, "FIXED_IRRIGATION");

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "DOY", &q->opDay, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "VOLUME", &q->opVolume, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                q->status = 0;
            }
        }

        cropmgmt->numAutoIrrigation = nautoirrg;
        if (nautoirrg > 0)
        {
            cropmgmt->autoIrrigation =
                (autoirr_struct *)malloc (nautoirrg *
                sizeof (autoirr_struct));
            /* Rewind to the beginning of file and read all planting operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nautoirrg; j++)
            {
                FindLine (op_file, "AUTO_IRRIGATION");

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "CROP",
                        cropmgmt->autoIrrigation[j].cropName, 's'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "START_DAY",
                        &cropmgmt->autoIrrigation[j].startDay, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "STOP_DAY",
                        &cropmgmt->autoIrrigation[j].stopDay, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "WATER_DEPLETION",
                        &cropmgmt->autoIrrigation[j].waterDepletion, 'd'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }

                NextLine (op_file, cmdstr);
                if (!ReadKeyword (cmdstr, "LAST_SOIL_LAYER",
                        &cropmgmt->autoIrrigation[j].lastSoilLayer, 'i'))
                {
                    fprintf (stderr, "Error opening %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
            }
        }

        /* Link plating order and auto irrigation */
        for (j = 0; j < cropmgmt->totalCropsPerRotation; j++)
        {
            if (cropmgmt->usingAutoIrr &&
                cropmgmt->plantingOrder[j].usesAutoIrrigation == 1)
            {
                for (k = 0; k < nautoirrg; k++)
                {
                    if (strcmp (cropmgmt->plantingOrder[j].cropName,
                            cropmgmt->autoIrrigation[k].cropName) == 0)
                    {
                        cropmgmt->plantingOrder[j].usesAutoIrrigation = k;
                        break;
                    }
                }
                if (k >= nautoirrg)
                {
                    printf
                        ("ERROR: Cannot find the description of auto irrigation for %s!\n",
                        cropmgmt->plantingOrder[j].cropName);
                    exit (1);
                }
            }
            else
                cropmgmt->plantingOrder[j].usesAutoIrrigation = -1;
        }

        fclose (op_file);
    }
}

int CropExist (char *cropName, const croptbl_struct *croptbl)
{
    int             i;
    int             exist = 0;

    for (i = 0; i < croptbl->number; i++)
    {
        if (strcmp (cropName, croptbl->cropName[i]) == 0)
        {
            exist = 1;
            break;
        }
    }

    return (exist);
}
