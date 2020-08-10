using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AppleTurnover
{
    public class PeptideTurnoverObject
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; set; }
        public string Protein { get; set; }
        public string Proteoform { get; set; }
        public double[] Timepoints { get; set; }
        public double[] RelativeFractions { get; set; }
        public string[] Filenames { get; set; }
        public double[] Intensities { get; set; }
        public double Kbi { get; set; }
        public double Error { get; set; }
        public double TemporaryError { get; set; }
        public double TotalIntensity { get; set; }
        public double LowKbi { get; set; }
        public double HighKbi { get; set; }
        public int StartResidue { get; set; }
        public int EndResidue { get; set; }
        public Dictionary<int, string> ModDictionary { get; set; }
        public string FileName { get; set; }
        public double[] MonteCarloKbis { get; set; }

        //properties used to display data
        public string ErrorString { get { return Math.Round(Error, 6).ToString(); } }
        public double Halflife { get { return Math.Round(Math.Log(2) / Kbi, 1); } }
        public double CI { get { return Math.Round((Math.Log(2) / LowKbi) - (Math.Log(2) / HighKbi), 1); } }
        public string DisplayProteinOrProteoform { get; private set; } //used to dynamically display either the protein level or proteoform level information
        public string DisplayPeptideSequence { get; private set; }//used to dynamically display either the full or base peptide sequence

        public PeptideTurnoverObject(string fullSequence, double[] timepoints, double[] rFValues, string[] filenames, double[] intensities, double totalIntensity, string fileName, string protein = "", string proteoform = null)
        {
            if (fullSequence.Contains('"'))
            { fullSequence = fullSequence.Replace("\"", ""); } //weird bug that sometimes occurs in file loading
            FullSequence = fullSequence;
            DisplayPeptideSequence = fullSequence;
            BaseSequence = CleanSeq(fullSequence);

            Protein = protein;
            DisplayProteinOrProteoform = protein;
            if(proteoform!=null && proteoform.Contains("@"))
            { }
            Proteoform = proteoform ?? protein;

            //sort by timepoints
            for (int i = 1; i < timepoints.Length; i++)
            {
                int j = i;
                while (j > 0 && timepoints[j] < timepoints[j - 1])
                {
                    double temp = timepoints[j - 1];
                    timepoints[j - 1] = timepoints[j];
                    timepoints[j] = temp;

                    temp = rFValues[j - 1];
                    rFValues[j - 1] = rFValues[j];
                    rFValues[j] = temp;

                    string tempString = filenames[j - 1];
                    filenames[j - 1] = filenames[j];
                    filenames[j] = tempString;

                    double tempIntensity = intensities[j - 1];
                    intensities[j - 1] = intensities[j];
                    intensities[j] = tempIntensity;
                    j--;
                }
            }

            //Array.Sort(timepoints, rFValues);
            Timepoints = timepoints;
            RelativeFractions = rFValues;
            Filenames = filenames;
            Intensities = intensities;
            Error = double.PositiveInfinity;
            TemporaryError = Error;
            TotalIntensity = totalIntensity;
            Kbi = 0.4; //arbitrary starting value
            LowKbi = 0;
            HighKbi = 2;
            ModDictionary = new Dictionary<int, string>();
            FileName = fileName;
        }

        public PeptideTurnoverObject Copy(string proteoform)
        {
            return new PeptideTurnoverObject(FullSequence, Timepoints, RelativeFractions, Filenames, Intensities, TotalIntensity, FileName, Protein, proteoform);
        }

        public void UpdateError()
        {
            Error = TemporaryError;
        }

        //removes all modifications (in brackets) and returns the base sequence
        public static string CleanSeq(string seq)
        {
            string cleanedSeq = "";
            char[] seqArray = seq.ToCharArray();
            int numBrackets = 0; //there can be nested brackets, such as [Fe[III]]
            foreach (char amino_acid in seqArray)
            {
                if (amino_acid == ']')
                {
                    numBrackets--;
                }
                else if (amino_acid == '[')
                {
                    numBrackets++;
                }
                else if (numBrackets == 0)
                {
                    cleanedSeq += amino_acid;
                }
            }
            return cleanedSeq;
        }

        public void UpdateDisplayProteinOrProteoform(bool displayProtein, bool displayFullSequence)
        {
            DisplayProteinOrProteoform = displayProtein ? Protein : Proteoform;
            DisplayPeptideSequence = displayFullSequence ? FullSequence : BaseSequence;
        }
    }
}