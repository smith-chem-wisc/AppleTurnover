using System;
using System.Collections.Generic;
using System.Text;

namespace TurnoverGUI
{
    public class Settings
    {
        public int MinValidValuesTotal { get; set; }
        public int MinValidValuesPerTimepoint { get; set; }
        public bool UseBadRatios { get; set; }
        public SearchEngine UpstreamProgram {get;set;}
        public Settings(int minValidValuesTotal = 6, int minValidValuesPerTimepoint= 3, bool useBadRatios = false, SearchEngine engine = SearchEngine.MaxQuant)
        {
            MinValidValuesTotal = minValidValuesTotal;
            MinValidValuesPerTimepoint = minValidValuesPerTimepoint;
            UseBadRatios = useBadRatios;
            UpstreamProgram = engine;
        }

        public enum SearchEngine
        {
            MaxQuant,
            MetaMorpheus
        }
    }
}
