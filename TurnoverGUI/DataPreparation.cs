﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Linq;
using Proteomics;
using UsefulProteomicsDatabases;
using System.Diagnostics;
using System.Threading.Tasks;
using FlashLFQ;
using System.Collections.ObjectModel;
using MathNet.Numerics.Statistics;
using Accessibility;
using ThermoFisher.CommonCore.Data;

namespace AppleTurnover
{
    public static class DataPreparation
    {
        public static readonly string UnmodifiedString = "Unmodified";
        public static List<PeptideTurnoverObject> ReadData(string file, Settings settings, List<Protein> theoreticalProteins)
        {
            List<string> columnIndexToSampleName = new List<string>();
            List<double> columnIndexToTimepoint = new List<double>();
            Dictionary<string, double> sampleToTimepoint = new Dictionary<string, double>();
            Dictionary<double, List<string>> timepointToSamples = new Dictionary<double, List<string>>();

            List<PeptideTurnoverObject> peptides = new List<PeptideTurnoverObject>();

            if (settings.UpstreamProgram == Settings.SearchEngine.MetaMorpheus)
            {
                int firstIntensityIndex = -1;
                string[] lines = File.ReadAllLines(file);
                string[] header = lines[0].Split('\t');
                //Get the header information for which columns are which samples and which timepoints
                for (int i = 0; i < header.Length; i += 2)
                {
                    if (!header[i].Contains("Intensity_"))
                    {
                        i++;
                        while (i < header.Length && !header[i].Contains("Intensity_"))
                        {
                            i++;
                        }
                        if (firstIntensityIndex == -1)
                        {
                            firstIntensityIndex = i;
                        }
                        i -= 2;
                    }
                    else
                    {
                        string columnHeader = header[i].Replace('-', '_');
                        string[] columnHeaderSplit = columnHeader.Split('_');
                        string sampleName = columnHeader.Substring(0, columnHeader.Length - columnHeaderSplit[columnHeaderSplit.Length - 1].Length - 1);
                        columnIndexToSampleName.Add(sampleName);
                        for (int j = 0; j < columnHeaderSplit.Length; j++)
                        {
                            string timepointString = columnHeaderSplit[j];
                            double time = -1;
                            if (timepointString[0] == 'd' || timepointString[0] == 'D')
                            {
                                try
                                {
                                    time = Convert.ToDouble(timepointString.Substring(1));
                                }
                                catch { };
                            }
                            else if (timepointString[timepointString.Length - 1] == 'd' || timepointString[timepointString.Length - 1] == 'D')
                            {
                                try
                                {
                                    time = Convert.ToDouble(timepointString.Substring(0, timepointString.Length - 1));
                                }
                                catch { };
                            }
                            if (time != -1)
                            {
                                if (timepointToSamples.ContainsKey(time))
                                {
                                    timepointToSamples[time].Add(sampleName);
                                }
                                else
                                {
                                    timepointToSamples[time] = new List<string> { sampleName };
                                }
                                sampleToTimepoint[sampleName] = time;
                                columnIndexToTimepoint.Add(time);
                                break;
                            }
                        }
                    }
                }

                //Diagnostics
                //int numOriginalPeps = 0;
                //int numSurvivingPeps = 0;
                ////List<string> originalProteins = new List<string>();
                ////List<string> survivingProteins = new List<string>();
                //int[] originalDistribution = new int[25];
                //int[] postBadRemovalDistribution = new int[25];
                //int[] postRequirementPerTimepoint = new int[25];
                //HashSet<string> originalProteins = new HashSet<string>();
                //HashSet<string> survivingProteins = new HashSet<string>();

                //10 investigate
                //int[] num3 = new int[6];
                //int[] num7 = new int[6];
                //int[] num14 = new int[6];
                //int[] num30 = new int[6];
                //int[] num60 = new int[6];
                //List<(double, string)> peptidesWith10Values = new List<(double, string)>();


                //read in the intensities
                bool peptideInput = header[0].Equals("Sequence");
                for (int i = 1; i < lines.Length; i++)
                {
                    string[] line = lines[i].Split('\t');
                    string sequence = line[0];
                    string protein = peptideInput ? line[2] : line[0];
                    List<double> timepointsForThisPeptide = new List<double>();
                    List<double> rfValuesForThisPeptide = new List<double>();
                    List<string> filenamesForThisPeptide = new List<string>();
                    List<double> intensitiesForThisPeptide = new List<double>();

                    //Diagnostics
                    //numOriginalPeps++;
                    //originalProteins.Add(protein);
                    int numOriginal = 0;
                    int numPostBadRemoval = 0;
                    int numPostRequirementPerTimepoint = 0;
                    double highestIntensity = 0;
                    double averageIntensity = 0;

                    for (int column = firstIntensityIndex; column < firstIntensityIndex + columnIndexToSampleName.Count * 2; column += 2)
                    {
                        double originalIntensity = line[column].Length == 0 ? 0 : Convert.ToDouble(line[column]);
                        double newlySynthesizedIntensity = line[column + 1].Length == 0 ? 0 : Convert.ToDouble(line[column + 1]);
                        bool atLeastOneIntensity = originalIntensity != 0 || newlySynthesizedIntensity != 0;
                        if (atLeastOneIntensity)
                        {
                            numOriginal++;
                            averageIntensity += originalIntensity + newlySynthesizedIntensity;
                            if (originalIntensity+newlySynthesizedIntensity>highestIntensity)
                            {
                                highestIntensity = originalIntensity + newlySynthesizedIntensity;
                            }
                        }
                        if ((settings.UseBadRatios && atLeastOneIntensity) || (originalIntensity != 0 && newlySynthesizedIntensity != 0))
                        {
                            numPostBadRemoval++;
                            int indexLookup = (column - firstIntensityIndex) / 2;
                            timepointsForThisPeptide.Add(columnIndexToTimepoint[indexLookup]);
                            rfValuesForThisPeptide.Add(newlySynthesizedIntensity / (originalIntensity + newlySynthesizedIntensity));
                            filenamesForThisPeptide.Add(columnIndexToSampleName[indexLookup]);
                            intensitiesForThisPeptide.Add(originalIntensity + newlySynthesizedIntensity);
                        }
                    }
                    //if (numOriginal == 10)
                    //{
                    //    num3[timepointsForThisPeptide.Count(x => x > 2 && x < 4)]++;
                    //    num7[timepointsForThisPeptide.Count(x => x > 6 && x < 8)]++;
                    //    num14[timepointsForThisPeptide.Count(x => x > 13 && x < 15)]++;
                    //    num30[timepointsForThisPeptide.Count(x => x > 29 && x < 31)]++;
                    //    num60[timepointsForThisPeptide.Count(x => x > 59 && x < 61)]++;
                    //    peptidesWith10Values.Add((intensitiesForThisPeptide.Median(), sequence));
                    //}
                    //remove timepoints with too little data
                    foreach (double timepoint in timepointToSamples.Keys)
                    {
                        List<int> indicesForThisTimepoint = new List<int>();
                        for (int index = 0; index < timepointsForThisPeptide.Count; index++)
                        {
                            if (timepointsForThisPeptide[index].Equals(timepoint))
                            {
                                numPostRequirementPerTimepoint++;
                                indicesForThisTimepoint.Add(index);
                            }
                        }
                        if (indicesForThisTimepoint.Count < settings.MinValidValuesPerTimepoint)
                        {
                            for (int index = indicesForThisTimepoint.Count - 1; index >= 0; index--)
                            {
                                numPostRequirementPerTimepoint--;
                                int actualIndex = indicesForThisTimepoint[index];
                                timepointsForThisPeptide.RemoveAt(actualIndex);
                                rfValuesForThisPeptide.RemoveAt(actualIndex);
                                filenamesForThisPeptide.RemoveAt(actualIndex);
                                intensitiesForThisPeptide.RemoveAt(actualIndex);
                            }
                        }
                    }

                    if (timepointsForThisPeptide.Count >= settings.MinValidValuesTotal)
                    {
                    //    numSurvivingPeps++;
                    //    survivingProteins.Add(protein);
                        peptides.Add(new PeptideTurnoverObject(sequence, timepointsForThisPeptide.ToArray(), rfValuesForThisPeptide.ToArray(),
                            filenamesForThisPeptide.ToArray(), intensitiesForThisPeptide.ToArray(), intensitiesForThisPeptide.Sum(), file, protein));
                    }
                    //originalDistribution[numOriginal]++;

                    //postBadRemovalDistribution[numPostBadRemoval]++;
                    //postRequirementPerTimepoint[numPostRequirementPerTimepoint]++;

                }
                string pathToWrite = file.Substring(0, file.Length - 4) + "_Results";
                Directory.CreateDirectory(pathToWrite);

                //DIAGNOSTICS
                //List<string> linesToWrite = new List<string>();
                //linesToWrite.Add("Original Peptides:\t" + numOriginalPeps.ToString());
                //linesToWrite.Add("Surviving Peptides:\t" + numSurvivingPeps.ToString());
                //linesToWrite.Add("Original Proteins:\t" + originalProteins.Count.ToString());
                //linesToWrite.Add("Surviving Proteins:\t" + survivingProteins.Count.ToString());
                //linesToWrite.Add("");
                //linesToWrite.Add("NumPeptides\tOriginal\tPostBadRatios\tPostTimepointMin");
                //for (int i = 0; i < 25; i++)
                //{
                //    linesToWrite.Add(i.ToString() + '\t' + originalDistribution[i].ToString() + '\t' + postBadRemovalDistribution[i].ToString() + '\t' + postRequirementPerTimepoint[i].ToString());
                //}
                //string filename = Path.GetFileNameWithoutExtension(file);

                //pathToWrite = Path.Combine(pathToWrite, filename);
                //File.WriteAllLines(pathToWrite + "_Diagnostics.tsv", linesToWrite);

                //linesToWrite.Clear();
                //linesToWrite.Add("\t3\t7\t14\t30\t60");
                //for(int i=0; i<6; i++)
                //{
                //    linesToWrite.Add(i.ToString() + '\t' + num3[i].ToString() + '\t' + num7[i].ToString() + '\t' + num14[i].ToString() + '\t' + num30[i].ToString() + '\t' + num60[i].ToString());
                //}
                //peptidesWith10Values = peptidesWith10Values.OrderByDescending(x => x.Item1).ToList();
                //linesToWrite.Add("");
                //linesToWrite.Add(peptidesWith10Values[0].Item2);
                //linesToWrite.Add(peptidesWith10Values[1].Item2);
                //linesToWrite.Add(peptidesWith10Values[2].Item2);
                //linesToWrite.Add(peptidesWith10Values[3].Item2);
                //linesToWrite.Add(peptidesWith10Values[4].Item2);
                //File.WriteAllLines(pathToWrite + "_DiagnosticsFor10ValidFiles.tsv", linesToWrite);
                //int a = 0;
            }
            else //if maxquant (not really maxquant, but the same format used in Alevra, M.; Mandad, S.; Ischebeck, T.; Urlaub, H.; Rizzoli, S. O.; Fornasiero, E. F. A Mass Spectrometry Workflow for Measuring Protein Turnover Rates in Vivo. Nature Protocols 2019. https://doi.org/10.1038/s41596-019-0222-y.)
            {
                int firstRatioIndex = -1;
                int firstIntensityIndex = -1;
                string[] lines = File.ReadAllLines(file);
                string[] header = lines[0].Split('\t');

                bool peptideInput = header[0].Equals("Sequence");

                //Get the header information for which columns are which samples and which timepoints
                if (peptideInput)
                {
                    for (int i = 0; i < header.Length; i += 6) //6 is the spacing for the actual, normalized, variablitiy, count, iso, type
                    {
                        if (!header[i].Contains("Ratio H/L"))
                        {
                            i++;
                            while (i < header.Length && !header[i].Contains("Ratio H/L"))
                            {
                                if (header[i].Contains("Intensity")) //intensities are after the ratios. Reason unknown, they do not match the ratios
                                {
                                    firstIntensityIndex = i + 3;
                                    i = header.Length; //end
                                }
                                i++;
                            }
                            if (firstRatioIndex == -1)
                            {
                                firstRatioIndex = i + 6;
                            }
                            //skip the first one because it's the aggregate
                        }
                        else
                        {
                            string sampleName = header[i].Replace('-', '_').Replace(' ', '_').Substring("Ratio H/L ".Length);
                            string[] columnHeaderSplit = sampleName.Split('_');
                            columnIndexToSampleName.Add(sampleName);

                            double time = Convert.ToDouble(columnHeaderSplit[0]);
                            if (timepointToSamples.ContainsKey(time))
                            {
                                timepointToSamples[time].Add(sampleName);
                            }
                            else
                            {
                                timepointToSamples[time] = new List<string> { sampleName };
                            }
                            sampleToTimepoint[sampleName] = time;
                            columnIndexToTimepoint.Add(time);
                        }
                    }

                    //read in the intensities
                    for (int i = 1; i < lines.Length; i++)
                    {
                        string[] line = lines[i].Split('\t');
                        string sequence = line[0];
                        string protein = line[34];
                        List<double> timepointsForThisPeptide = new List<double>();
                        List<double> rfValuesForThisPeptide = new List<double>();
                        List<string> filenamesForThisPeptide = new List<string>();
                        List<double> intensitiesForThisPeptide = new List<double>();
                        for (int column = firstRatioIndex; column < firstRatioIndex + columnIndexToSampleName.Count * 6; column += 6)
                        {
                            if (!line[column].Equals("NaN"))
                            {
                                int indexLookup = (column - firstRatioIndex) / 6;
                                timepointsForThisPeptide.Add(columnIndexToTimepoint[indexLookup]);

                                double ratio = Convert.ToDouble(line[column]);
                                rfValuesForThisPeptide.Add(1 - ratio / (ratio + 1)); //convert H/L to L/Total //TODO remove the "1-" for normal experiments before release
                                filenamesForThisPeptide.Add(columnIndexToSampleName[indexLookup]);
                                intensitiesForThisPeptide.Add(Convert.ToDouble(line[(column - firstRatioIndex) / 6 + firstIntensityIndex]));
                            }
                        }
                        //remove timepoints with too little data
                        foreach (double timepoint in timepointToSamples.Keys)
                        {
                            List<int> indicesForThisTimepoint = new List<int>();
                            for (int index = 0; index < timepointsForThisPeptide.Count; index++)
                            {
                                if (timepointsForThisPeptide[index].Equals(timepoint))
                                {
                                    indicesForThisTimepoint.Add(index);
                                }
                            }
                            if (indicesForThisTimepoint.Count < settings.MinValidValuesPerTimepoint)
                            {
                                for (int index = indicesForThisTimepoint.Count - 1; index >= 0; index--)
                                {
                                    int actualIndex = indicesForThisTimepoint[index];
                                    timepointsForThisPeptide.RemoveAt(actualIndex);
                                    rfValuesForThisPeptide.RemoveAt(actualIndex);
                                    filenamesForThisPeptide.RemoveAt(actualIndex);
                                    intensitiesForThisPeptide.RemoveAt(actualIndex);
                                }
                            }
                        }

                        if (timepointsForThisPeptide.Count >= settings.MinValidValuesTotal)
                        {
                            peptides.Add(new PeptideTurnoverObject(sequence, timepointsForThisPeptide.ToArray(), rfValuesForThisPeptide.ToArray(),
                                filenamesForThisPeptide.ToArray(), intensitiesForThisPeptide.ToArray(), intensitiesForThisPeptide.Sum(), file, protein));
                        }
                    }
                }
                else //if protein input
                {

                    for (int i = 0; i < header.Length; i++) //no spacing here
                    {
                        if (!header[i].Contains("Ratio H/L"))
                        {
                            i++;
                            while (i < header.Length && !header[i].Contains("Ratio H/L"))
                            {
                                if (header[i].Contains("Intensity")) //intensities are after the ratios. Reason unknown, they do not match the ratios
                                {
                                    firstIntensityIndex = i + 3;
                                    i = header.Length; //end
                                }
                                i++;
                            }
                            if (firstRatioIndex == -1)
                            {
                                firstRatioIndex = i;
                            }
                            //skip the first one because it's the aggregate
                        }
                        else
                        {
                            if (header[i].Contains("normalized"))
                            {
                                while (i < header.Length)
                                {
                                    if (header[i].Contains("Intensity")) //intensities are after the ratios. Reason unknown, they do not match the ratios
                                    {
                                        firstIntensityIndex = i + 3;
                                        i = header.Length; //end
                                    }
                                    i++;
                                }
                                break;
                            }
                            string sampleName = header[i].Replace('-', '_').Replace(' ', '_').Substring("Ratio H/L ".Length);
                            string[] columnHeaderSplit = sampleName.Split('_');
                            columnIndexToSampleName.Add(sampleName);

                            double time = Convert.ToDouble(columnHeaderSplit[0]);
                            if (timepointToSamples.ContainsKey(time))
                            {
                                timepointToSamples[time].Add(sampleName);
                            }
                            else
                            {
                                timepointToSamples[time] = new List<string> { sampleName };
                            }
                            sampleToTimepoint[sampleName] = time;
                            columnIndexToTimepoint.Add(time);
                        }
                    }

                    //read in the intensities
                    for (int i = 1; i < lines.Length; i++)
                    {
                        string[] line = lines[i].Split('\t');
                        string sequence = line[0];
                        string protein = line[0];
                        List<double> timepointsForThisPeptide = new List<double>();
                        List<double> rfValuesForThisPeptide = new List<double>();
                        List<string> filenamesForThisPeptide = new List<string>();
                        List<double> intensitiesForThisPeptide = new List<double>();
                        for (int column = firstRatioIndex; column < firstRatioIndex + columnIndexToSampleName.Count; column++)
                        {
                            if (!line[column].Equals("NaN"))
                            {
                                int indexLookup = (column - firstRatioIndex);
                                timepointsForThisPeptide.Add(columnIndexToTimepoint[indexLookup]);

                                double ratio = Convert.ToDouble(line[column]);
                                rfValuesForThisPeptide.Add(1 - ratio / (ratio + 1)); //convert H/L to L/Total //TODO remove the "1-" for normal experiments before release
                                filenamesForThisPeptide.Add(columnIndexToSampleName[indexLookup]);
                                intensitiesForThisPeptide.Add(Convert.ToDouble(line[(column - firstRatioIndex) + firstIntensityIndex]));
                            }
                        }
                        //remove timepoints with too little data
                        foreach (double timepoint in timepointToSamples.Keys)
                        {
                            List<int> indicesForThisTimepoint = new List<int>();
                            for (int index = 0; index < timepointsForThisPeptide.Count; index++)
                            {
                                if (timepointsForThisPeptide[index].Equals(timepoint))
                                {
                                    indicesForThisTimepoint.Add(index);
                                }
                            }
                            if (indicesForThisTimepoint.Count < settings.MinValidValuesPerTimepoint)
                            {
                                for (int index = indicesForThisTimepoint.Count - 1; index >= 0; index--)
                                {
                                    int actualIndex = indicesForThisTimepoint[index];
                                    timepointsForThisPeptide.RemoveAt(actualIndex);
                                    rfValuesForThisPeptide.RemoveAt(actualIndex);
                                    filenamesForThisPeptide.RemoveAt(actualIndex);
                                    intensitiesForThisPeptide.RemoveAt(actualIndex);
                                }
                            }
                        }

                        if (timepointsForThisPeptide.Count >= settings.MinValidValuesTotal)
                        {
                            peptides.Add(new PeptideTurnoverObject(sequence, timepointsForThisPeptide.ToArray(), rfValuesForThisPeptide.ToArray(),
                                filenamesForThisPeptide.ToArray(), intensitiesForThisPeptide.ToArray(), intensitiesForThisPeptide.Sum(), file, protein));
                        }
                    }
                }
            }

            //create a hash table for quick lookups
            //the idea is to break up the proteins into k-mers (must be less than shortest peptide length) and look up the starts of those

            Dictionary<Protein, Dictionary<string, List<int>>> theoreticalProteinLookupTable = new Dictionary<Protein, Dictionary<string, List<int>>>();
            const int kMerLength = 6;
            foreach (Protein p in theoreticalProteins)
            {
                string baseSequence = p.BaseSequence;
                Dictionary<string, List<int>> lookupTable = new Dictionary<string, List<int>>();
                for (int i = 0; i < baseSequence.Length - kMerLength + 1; i++)
                {
                    string kmer = baseSequence.Substring(i, kMerLength);
                    if (lookupTable.TryGetValue(kmer, out var value))
                    {
                        value.Add(i);
                    }
                    else
                    {
                        lookupTable.Add(kmer, new List<int> { i });
                    }
                }
                if (!theoreticalProteinLookupTable.Keys.Contains(p))
                {
                    theoreticalProteinLookupTable.Add(p, lookupTable);
                }
            }

            //lookup the sequence in the database
            //sort
            peptides = peptides.Where(x => x.RelativeFractions.Length != 0).OrderBy(x => x.BaseSequence).ToList();

            //do reverse parsimony
            int[] threads = Enumerable.Range(0, Environment.ProcessorCount).ToArray();
            //int[] threads = Enumerable.Range(0, 1).ToArray();
            Parallel.ForEach(threads, (thread) =>
            {
                string mostRecentBaseSequence = "";
                int mostRecentIndex = -1;
                int max = (thread + 1) * peptides.Count / threads.Length;
                for (int i = thread * peptides.Count / threads.Length; i < max; i++)
                {
                    PeptideTurnoverObject currentPeptide = peptides[i];
                    //if same base seq, just reuse the old info
                    if (!currentPeptide.BaseSequence.Equals(mostRecentBaseSequence))
                    {
                        mostRecentBaseSequence = currentPeptide.BaseSequence;

                        //find proteins containing this sequence
                        //List<Protein> proteinsContainingThisSeq = theoreticalProteins.Where(x => x.BaseSequence.Contains(mostRecentBaseSequence)).OrderBy(x => x.Accession).ToList();
                        List<Protein> proteinsContainingThisSeq = new List<Protein>();

                        string kMer = mostRecentBaseSequence.Substring(0, kMerLength);
                        foreach (Protein p in theoreticalProteins)
                        {
                            if (theoreticalProteinLookupTable[p].TryGetValue(kMer, out var indices))
                            {
                                if (NativeStringSearch(p.BaseSequence, mostRecentBaseSequence, indices, kMer))
                                {
                                    proteinsContainingThisSeq.Add(p);
                                }
                            }
                        }
                        if (proteinsContainingThisSeq.Count == 0)
                        {
                            throw new Exception("Database was missing the protein: " + currentPeptide.Protein + " or the given protein did not contain the sequence: " + mostRecentBaseSequence);
                        }
                        string protein = proteinsContainingThisSeq[0].Accession;
                        for (int index = 1; index < proteinsContainingThisSeq.Count; index++)
                        {
                            protein += ";" + proteinsContainingThisSeq[index].Accession;
                        }
                        if (!protein.Equals(currentPeptide.Protein))
                        { //why aren't they the same? Parsimony
                            currentPeptide.UpdateProteinFromParsimony(protein);
                        }

                        if (proteinsContainingThisSeq.Count > 1)
                        { } //TODO: The handling of this is sloppy, but complicated to deal with the right way. 
                            //It's assuming that proteoforms will only appear by being on the same peptide (modified/unmodified) and not through overlapping peptides
                            //find index of base sequence 
                        mostRecentIndex = proteinsContainingThisSeq.First().BaseSequence.IndexOf(mostRecentBaseSequence);
                    }
                    currentPeptide.StartResidue = mostRecentIndex;
                    currentPeptide.EndResidue = mostRecentIndex + mostRecentBaseSequence.Length;
                    //if there are mods
                    if (mostRecentBaseSequence.Length != currentPeptide.FullSequence.Length)
                    {
                        string fullSequence = currentPeptide.FullSequence;
                        int currentIndex = 0;
                        for (int index = 0; index < fullSequence.Length; index++)
                        {
                            //if there's a mod
                            if (fullSequence[index] == '[')
                            {
                                int bracketCount = 1;
                                string mod = "";
                                index++;
                                while (bracketCount != 0)
                                {
                                    if (fullSequence[index] == '[')
                                    {
                                        bracketCount++;
                                    }
                                    else if (fullSequence[index] == ']')
                                    {
                                        bracketCount--;
                                        index--;
                                    }
                                    if (bracketCount != 0)
                                    {
                                        mod += fullSequence[index];
                                    }
                                    index++;
                                }
                                if (currentIndex == 1 && currentPeptide.ModDictionary.ContainsKey(0))
                                {
                                    currentPeptide.ModDictionary[0] = currentPeptide.ModDictionary[0] + "+" + mod;
                                }
                                else
                                {
                                    currentPeptide.ModDictionary.Add(currentIndex == 0 ? 0 : currentIndex - 1, mod); //N-terminal mods are counted as being on the first residue
                                }
                            }
                            else
                            {
                                currentIndex++;
                            }
                        }
                    }
                }
            });

            //have all peptides, now convert into proteins
            List<PeptideTurnoverObject> possibleProteoformGroups = new List<PeptideTurnoverObject>();
            var peptidesGroupedByProtein = peptides.GroupBy(x => x.Protein, x => x).ToList();

            //TODO parallelize this, break out into separate method, unit tests
            //foreach protein group
            for (int i = 0; i < peptidesGroupedByProtein.Count; i++)
            {
                var group = peptidesGroupedByProtein[i];
                string protein = group.Key;
                //get the peptides
                List<PeptideTurnoverObject> peptidesForThisProtein = group.OrderBy(x => x.StartResidue).ToList();
                Dictionary<int, List<string>> mods = new Dictionary<int, List<string>>();
                foreach (PeptideTurnoverObject peptide in peptidesForThisProtein)
                {
                    for (int residue = peptide.StartResidue; residue < peptide.EndResidue; residue++)
                    {
                        //is residue start/end correct?
                        string value = UnmodifiedString;
                        if (peptide.ModDictionary.ContainsKey(residue - peptide.StartResidue)) //the mod dictionary is zero-based.
                        {
                            value = peptide.ModDictionary[residue - peptide.StartResidue];
                        }
                        if (mods.ContainsKey(residue))
                        {
                            if (!mods[residue].Contains(value))
                            {
                                mods[residue].Add(value);
                            }
                        }
                        else
                        {
                            mods.Add(residue, new List<string> { value });
                        }
                    }
                }

                //find which residues have multiple forms
                var residueDifferences = mods.Where(x => x.Value.Count > 1).OrderBy(x => x.Key).ToList();
                bool[] uniquePeptides = new bool[peptidesForThisProtein.Count]; //have we added this peptide already? At the start, none have been added.
                if (residueDifferences.Count != 0)
                {
                    //foreach residue with a different form
                    foreach (var residueDifference in residueDifferences)
                    {
                        int residueIndex = residueDifference.Key; //get the residue

                        //find the relevant peptides for this residue
                        bool foundResidue = false;
                        for (int index = 0; index < peptidesForThisProtein.Count; index++)
                        {
                            PeptideTurnoverObject peptide = peptidesForThisProtein[index];
                            //if this peptide is relevant as belonging to a proteoform group
                            if (peptide.StartResidue <= residueIndex && peptide.EndResidue >= residueIndex)
                            {//is this right?
                                uniquePeptides[index] = true;
                                foundResidue = true;
                                string modForThisPeptide = UnmodifiedString;
                                if (peptide.ModDictionary.ContainsKey(residueIndex - peptide.StartResidue))
                                {
                                    modForThisPeptide = peptide.ModDictionary[residueIndex - peptide.StartResidue];
                                }
                                for (int indexForThisMod = 0; indexForThisMod < residueDifference.Value.Count; indexForThisMod++)
                                {
                                    if (residueDifference.Value[indexForThisMod].Equals(modForThisPeptide))
                                    {
                                        possibleProteoformGroups.Add(peptide.Copy(peptide.Protein + "_" + modForThisPeptide + "@" + residueDifference.Key.ToString()));
                                        break;
                                    }
                                }
                            }
                            //if we're no longer looking at relevant peptides
                            else if (foundResidue)
                            {
                                break;
                            }
                        }
                    }
                }

                //add peptides that aren't part of proteoform groups
                for (int index = 0; index < peptidesForThisProtein.Count; index++)
                {
                    if (!uniquePeptides[index])
                    {
                        possibleProteoformGroups.Add(peptidesForThisProtein[index]);
                    }
                }
            }

            return possibleProteoformGroups; //.Where(x => x.Timepoints.Length >= settings.MinValidValuesTotal).OrderByDescending(x => x.Timepoints.Length).ThenByDescending(x => x.TotalIntensity).ToList();
        }

        //Remove peptides that are shared in multiple proteins
        public static void RemoveSharedPeptides(List<PeptideTurnoverObject> peptides, List<Protein> proteins)
        {
            int peptidesRemoved = 0;
            List<string> proteinSequences = proteins.Select(x => x.BaseSequence).ToList();
            for (int i = peptides.Count - 1; i > 0; i--)
            {
                string peptideBaseSequence = peptides[i].BaseSequence;
                if (proteinSequences.Count(x => x.Contains(peptideBaseSequence)) > 1)
                {
                    peptides.RemoveAt(i);
                    peptidesRemoved++;
                }
            }
        }

        public static List<Protein> LoadProteins(List<string> dbFilenameList)
        {
            List<Protein> proteinList = new List<Protein>();
            foreach (var db in dbFilenameList)
            {
                var dbProteinList = LoadProteinDb(db);
                proteinList = proteinList.Concat(dbProteinList).ToList();
            }
            return proteinList;
        }

        private static bool NativeStringSearch(string target, string query, List<int> indices, string kMer)
        {
            foreach (int index in indices)
            {
                if (index + query.Length <= target.Length)
                {
                    for (int i = kMer.Length; i < query.Length; i++)
                    {
                        if (target[index + i] == query[i])
                        {
                            if (i == query.Length - 1)
                            {
                                return true;
                            }
                        }
                        else { break; }
                    }
                }
            }
            return false;
        }

        private static List<Protein> LoadProteinDb(string fileName)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                proteinList = ProteinDbLoader.LoadProteinFasta(fileName, true, DecoyType.None, false, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, out dbErrors);
            }
            else
            {
                proteinList = ProteinDbLoader.LoadProteinXML(fileName, true, DecoyType.None, null, false, null, out var um);
            }

            return proteinList.Where(p => p.BaseSequence.Length > 0).ToList();
        }

        public static bool LoadExistingResults(string inputFile, string fileToLoad, Dictionary<string, PoolParameters> poolParameterDictionary, ObservableCollection<PeptideTurnoverObject> peptides, List<PeptideTurnoverObject> proteins, List<PeptideTurnoverObject> proteoforms)
        {
            try
            {
                string[] lines = File.ReadAllLines(fileToLoad);
                //We need to read in the:
                //-pool parameters
                double[] poolParams = lines[0].Split('\t').Select(x => Convert.ToDouble(x)).ToArray();
                poolParameterDictionary[inputFile] = new PoolParameters(poolParams[0], poolParams[1], poolParams[2]);

                //-peptides
                int i = 1;
                for (; i < lines.Length; i++)
                {
                    string[] line = lines[i].Split('\t').ToArray();
                    if (line.Length == 1)
                    {
                        i++;
                        break;
                    }
                    PeptideTurnoverObject peptide = new PeptideTurnoverObject(
                        line[0],
                        line[1].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        line[2].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        line[3].Split(';'),
                        line[4].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        Convert.ToDouble(line[5]),
                        inputFile, //file
                        line[7],
                        line[8]);
                    peptide.Kbi = Convert.ToDouble(line[9]);
                    peptide.Error = Convert.ToDouble(line[10]);
                    peptide.LowKbi = Convert.ToDouble(line[11]);
                    peptide.HighKbi = Convert.ToDouble(line[12]);
                    peptides.Add(peptide);

                }
                //-Proteins
                for (; i < lines.Length; i++)
                {
                    string[] line = lines[i].Split('\t').ToArray();
                    if (line.Length == 1)
                    {
                        i++;
                        break;
                    }
                    PeptideTurnoverObject protein = new PeptideTurnoverObject(
                        line[0],
                        line[1].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        line[2].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        line[3].Split(';'),
                        line[4].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        Convert.ToDouble(line[5]),
                        inputFile, //file
                        line[7]);
                    protein.Kbi = Convert.ToDouble(line[8]);
                    protein.LowKbi = Convert.ToDouble(line[9]);
                    protein.HighKbi = Convert.ToDouble(line[10]);
                    proteins.Add(protein);
                }
                //-Proteoforms
                for (; i < lines.Length; i++)
                {
                    string[] line = lines[i].Split('\t').ToArray();

                    PeptideTurnoverObject proteoform = new PeptideTurnoverObject(
                        line[0],
                        line[1].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        line[2].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        line[3].Split(';'),
                        line[4].Split(';').Select(x => Convert.ToDouble(x)).ToArray(),
                        Convert.ToDouble(line[5]),
                        inputFile, //file
                        line[7]);
                    proteoform.Kbi = Convert.ToDouble(line[8]);
                    proteoform.LowKbi = Convert.ToDouble(line[9]);
                    proteoform.HighKbi = Convert.ToDouble(line[10]);
                    proteoforms.Add(proteoform);
                }
                return true;
            }
            catch
            {
                return false;
            }
        }

        public static void WriteQuickLoadFile(string fileToWrite, PoolParameters poolParams, List<PeptideTurnoverObject> peptides, List<PeptideTurnoverObject> proteins, List<PeptideTurnoverObject> proteoforms)
        {
            List<string> linesToWrite = new List<string>();
            //We need to write the:
            //-pool parameters
            linesToWrite.Add(poolParams.Kst.ToString() + '\t' + poolParams.Kbt.ToString() + '\t' + poolParams.Kao.ToString());
            //-Peptides
            foreach (PeptideTurnoverObject peptide in peptides)
            {
                linesToWrite.Add(
                    peptide.FullSequence + '\t' +
                    string.Join(';', peptide.Timepoints) + '\t' +
                    string.Join(';', peptide.RelativeFractions) + '\t' +
                    string.Join(';', peptide.Filenames) + '\t' +
                    string.Join(';', peptide.Intensities) + '\t' +
                    peptide.TotalIntensity.ToString() + '\t' +
                    peptide.FileName + '\t' +
                    peptide.Protein + '\t' +
                    peptide.Proteoform + '\t' +
                    peptide.Kbi.ToString() + '\t' +
                    peptide.Error.ToString() + '\t' +
                    peptide.LowKbi.ToString() + '\t' +
                    peptide.HighKbi.ToString());
            }
            linesToWrite.Add("-");
            //-Proteins
            foreach (PeptideTurnoverObject protein in proteins)
            {
                linesToWrite.Add(
                    protein.FullSequence + '\t' +
                    string.Join(';', protein.Timepoints) + '\t' +
                    string.Join(';', protein.RelativeFractions) + '\t' +
                    string.Join(';', protein.Filenames) + '\t' +
                    string.Join(';', protein.Intensities) + '\t' +
                    protein.TotalIntensity.ToString() + '\t' +
                    protein.FileName + '\t' +
                    protein.Protein + '\t' +
                    protein.Kbi.ToString() + '\t' +
                    protein.LowKbi.ToString() + '\t' +
                    protein.HighKbi.ToString());
            }
            linesToWrite.Add("-");
            //-Proteoforms
            foreach (PeptideTurnoverObject proteoform in proteoforms)
            {
                linesToWrite.Add(
                    proteoform.FullSequence + '\t' +
                    string.Join(';', proteoform.Timepoints) + '\t' +
                    string.Join(';', proteoform.RelativeFractions) + '\t' +
                    string.Join(';', proteoform.Filenames) + '\t' +
                    string.Join(';', proteoform.Intensities) + '\t' +
                    proteoform.TotalIntensity.ToString() + '\t' +
                    proteoform.FileName + '\t' +
                    proteoform.Protein + '\t' +
                    proteoform.Kbi.ToString() + '\t' +
                    proteoform.LowKbi.ToString() + '\t' +
                    proteoform.HighKbi.ToString());
            }
            File.WriteAllLines(fileToWrite, linesToWrite);
        }
    }
}