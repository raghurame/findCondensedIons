
/* calculate average degree of chain ionization (non-exclusion method) from LAMMPS output */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

/*
 * ARGS TO PASS
 * ~~~~~~~~~~~~
 *
 * argv[0] = program
 * argv[1] = input file name (LAMMPS dump file with number, type, x, y, z as columns expected)
 * argv[2] = cutoff distance
 * argv[3] = timesteps to skip initially
 * argv[4] = Type of sulphur atom in polymer chain in LAMMPS input file
 * argv[5] = number of monomers in the polymer chain (assumes only one polymer chain)
 * argv[6] = Type of counterion in LAMMPS input file
 * argv[7] = number of counterions in the simulation
 * argv[8] = Valency of counterion
 */

// Use this function if the input coords are unscaled
float convertToScaled (float xlo, float xhi, float unScaledCoords)
{
	float scaledCoords;
	scaledCoords = (unScaledCoords - xlo) / (xhi - xlo);
	return scaledCoords;
}

typedef struct cartesianPosition
{
	double x, y, z;
} CARTESIAN;

typedef struct lammpsDump
{
	int sino, atomType, ix, iy, iz;
	double x, y, z;
} LAMMPSDUMP;

int main(int argc, char *argv[])
{
	int nMonomers = atoi (argv[5]), nCounterIons = atoi (argv[7]), nAtomsTotal = nMonomers + nCounterIons, *counterIonFlag, currentTimeframe = 0, currentLineNumber = 0, currentCounterIon = 0, currentBead = 0, beadAtomType = atoi (argv[4], counterIonAtomType = atoi (argv[6])), timeframesToSkip = atoi (argv[3]), nCondensedIons = 0, counterIonValency = atoi (argv[8]);
	counterIonFlag = (int *) calloc (nCounterIons, sizeof (int));

	double rCutoff = atoi (argv[2]), distance, alphaValue;

	char lineString[1000];

	FILE *inputFile, *outputFile, *outputLog;
	inputFile = fopen (argv[1], "r");
	outputFile = fopen ("nonexclusion.output", "w");
	outputLog = fopen ("nonexclusion.logs", "w");

	CARTESIAN *bead, *counterIon, position, positionDifference, lowerBounds, upperBounds;
	bead = (CARTESIAN *) malloc (nMonomers * sizeof (CARTESIAN));
	counterIon = (CARTESIAN *) malloc (nCounterIons * sizeof (CARTESIAN));

	LAMMPSDUMP inputDump;

	while ((fgets (lineString, 1000, inputFile) != NULL))
	{
		currentLineNumber++;

		if (currentLineNumber == 6)
			sscanf (lineString, "%lf %lf\n", lowerBounds.x, upperBounds.x);
		if (currentLineNumber == 7)
			sscanf (lineString, "%lf %lf\n", lowerBounds.y, upperBounds.y);
		if (currentLineNumber == 8)
			sscanf (lineString, "%lf %lf\n", lowerBounds.z, upperBounds.z);

		// First 9 lines in LAMMPS dump are headers
		if (currentLineNumber > 9 && currentLineNumber <= (9 + nAtomsTotal))
		{
			sscanf (lineString, "%d %d %lf %lf %lf", inputDump.sino, inputDump.atomType, inputDump.x, inputDump.y, inputDump.z);

			if (inputDump.atomType == beadAtomType)
			{
				// Input the info if necessary from the LAMMPS dump file
				inputDump.ix = 0;
				inputDump.iy = 0;
				inputDump.iz = 0;
				bead[currentBead].x = lowerBounds.x + (upperBounds.x - lowerBounds.x) * (inputDump.x + (double) inputDump.ix);
				bead[currentBead].y = lowerBounds.y + (upperBounds.y - lowerBounds.y) * (inputDump.y + (double) inputDump.iy);
				bead[currentBead].z = lowerBounds.z + (upperBounds.z - lowerBounds.z) * (inputDump.z + (double) inputDump.iz);
				currentBead++;
			}
			else if (inputDump.atomType == counterIonAtomType)
			{
				// Input the info if necessary from the LAMMPS dump file
				inputDump.ix = 0;
				inputDump.iy = 0;
				inputDump.iz = 0;
				counterIon[currentCounterIon].x = lowerBounds.x + (upperBounds.x - lowerBounds.x) * (inputDump.x + (double) inputDump.ix);
				counterIon[currentCounterIon].y = lowerBounds.y + (upperBounds.y - lowerBounds.y) * (inputDump.y + (double) inputDump.iy);
				counterIon[currentCounterIon].z = lowerBounds.z + (upperBounds.z - lowerBounds.z) * (inputDump.z + (double) inputDump.iz);
				currentCounterIon++;
			}
		}

		if (currentLineNumber == (9 + nAtomsTotal))
		{
			currentTimeframe++;

			// Reset counts to be used for next timeframe
			currentCounterIon = 0;
			currentBead = 0;
			currentLineNumber = 0;

			// Skip the first few timeframes pre-equilibration
			if (currentTimeframe >= timeframesToSkip)
			{
				for (int monomer = 0; i < nMonomers; ++i)
				{
					for (int ion = 0; i < nCounterIons; ++i)
					{
						positionDifference.x = bead[monomer].x - counterIon[ion].x;
						positionDifference.y = bead[monomer].y - counterIon[ion].y;
						positionDifference.z = bead[monomer].z - counterIon[ion].z;
						distance = sqrt (pow (positionDifference.x, 2) + pow (positionDifference.y, 2), pow (positionDifference.z, 2));

						if (distance <= rCutoff && counterIonFlag[ion] == 0)
							nCondensedIons++;
					}
				}
				fprintf (stdout, "Timeframe: %d; Number of condensed ions (per timeframe): %lf", currentTimeframe, (double) (nCondensedIons / currentTimeframe));
				fprintf (outputLog, "Timeframe: %d; Number of condensed ions (per timeframe): %lf", currentTimeframe, (double) (nCondensedIons / currentTimeframe));
				fflush (stdout);
			}
		}
	}

	totalTimeframes = currentTimeframe;
	alphaValue = (double) (nMonomers - (counterIonValency * nCondensedIons)) / totalTimeframes;
	alphaValue /= (double) nMonomers;
	fprintf (stdout, "alpha value: %lf", alphaValue);
	fprintf (outputFile, "alpha value: %lf", alphaValue);

	fclose (inputFile);
	fclose (outputFile);
	fclose (outputLog);

	return 0;
}









