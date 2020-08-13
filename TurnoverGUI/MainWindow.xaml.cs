using FlashLFQ;
using LiveCharts.Wpf;
using Proteomics;
using ScottPlot.Statistics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Threading;
using UsefulProteomicsDatabases;

namespace AppleTurnover
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RawDataForDataGrid> DataFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private readonly ObservableCollection<RawDataForDataGrid> DatabasesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        ICollectionView DisplayPeptidesView;
        private bool AllowFileDrop = true;
        private readonly ObservableCollection<PeptideTurnoverObject> PeptidesToDisplay = new ObservableCollection<PeptideTurnoverObject>();
        private readonly ObservableCollection<PeptideTurnoverObject> AllPeptides = new ObservableCollection<PeptideTurnoverObject>();
        private readonly ObservableCollection<string> FilesToDisplayObservableCollection = new ObservableCollection<string>();
        private readonly ObservableCollection<string> FilesToHideObservableCollection = new ObservableCollection<string>();
        private Dictionary<string, PoolParameters> PoolParameterDictionary = new Dictionary<string, PoolParameters>();
        private List<PeptideTurnoverObject> AnalyzedProteins = new List<PeptideTurnoverObject>();
        private List<PeptideTurnoverObject> AnalyzedProteoforms = new List<PeptideTurnoverObject>();
        private List<PeptideTurnoverObject> PeptidesPreviouslyPlotted = new List<PeptideTurnoverObject>();
        private bool ChangeParamTextBox;
        private static bool DisplayProteinInSpecificTable = true; //display the protein (true) or the proteoform (false)
        private static bool DisplayFullPeptideSequence;

        public MainWindow()
        {
            InitializeComponent();
            dataGridPeptideFiles.DataContext = DataFilesObservableCollection;
            dataGridDatabaseFiles.DataContext = DatabasesObservableCollection;
            DisplayedSamplesDataGrid.DataContext = FilesToDisplayObservableCollection;
            HiddenSamplesDataGrid.DataContext = FilesToHideObservableCollection;
            DisplayAnalyzedFilesDataGrid.DataContext = DataFilesObservableCollection;
            HalfLifeHistogramPlot.Configure(enableScrollWheelZoom: false);
            PrecisionPlot.Configure(enableScrollWheelZoom: false);
            HalfLifeComparisonPlot.Configure(enableScrollWheelZoom: false);
            RatioComparisonPlot.Configure(enableScrollWheelZoom: false);
            peptideRadioButton.IsChecked = true;
            proteinSpecificRadioButton.IsChecked = DisplayProteinInSpecificTable;

            DisplayPeptidesView = CollectionViewSource.GetDefaultView(PeptidesToDisplay);
            DisplayPeptidesDataGrid.DataContext = DisplayPeptidesView;
            Loaders.LoadElements();
            PopulateChoices();
        }

        private void PopulateChoices()
        {
            Settings defaultSettings = new Settings();
            MinValidValuesTotalTextBox.Text = defaultSettings.MinValidValuesTotal.ToString();
            MinValidValuesTimepointTextBox.Text = defaultSettings.MinValidValuesPerTimepoint.ToString();
            UseBadRatiosCheckBox.IsChecked = defaultSettings.UseBadRatios;
            MaxQuantRadioButton.IsChecked = defaultSettings.UpstreamProgram == Settings.SearchEngine.MaxQuant;
            MetaMorpheusRadioButton.IsChecked = defaultSettings.UpstreamProgram == Settings.SearchEngine.MetaMorpheus;
            RemoveBadPeptidesCheckBox.IsChecked = defaultSettings.RemoveMessyPeptides;
            showFullSequenceCheckbox.IsChecked = true;
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (AllowFileDrop)
            {
                string[] files = null;
                try
                {
                    files = ((string[])e.Data.GetData(DataFormats.FileDrop)).OrderBy(p => p).ToArray();
                }
                catch
                {
                    //do nothing; the user dragged something which is not a file (like text)
                }

                if (files != null)
                {
                    foreach (var draggedFilePath in files)
                    {
                        if (Directory.Exists(draggedFilePath))
                        {
                            foreach (string file in Directory.EnumerateFiles(draggedFilePath, "*.*", SearchOption.AllDirectories))
                            {
                                AddAFile(file);
                            }
                        }
                        else
                        {
                            AddAFile(draggedFilePath);
                        }
                    }

                    dataGridPeptideFiles.CommitEdit(DataGridEditingUnit.Row, true);
                    dataGridPeptideFiles.Items.Refresh();
                    dataGridDatabaseFiles.CommitEdit(DataGridEditingUnit.Row, true);
                    dataGridDatabaseFiles.Items.Refresh();
                }
            }
        }

        private void AddAFile(string filepath)
        {
            var theExtension = Path.GetExtension(filepath).ToLowerInvariant();

            switch (theExtension)
            {
                case ".tsv":
                case ".psmtsv":
                case ".txt":
                    filepath = filepath.Replace("_ApplETurnoverSavedSession.tsv", ".tsv"); //replace quickload, if applicable
                    if (!DataFilesObservableCollection.Any(x => x.FilePath.Equals(filepath)))
                    {
                        DataFilesObservableCollection.Add(new RawDataForDataGrid(filepath));
                    }

                    break;
                case ".fasta":
                case ".xml":
                    if (!DatabasesObservableCollection.Any(x => x.FilePath.Equals(filepath)))
                    {
                        DatabasesObservableCollection.Add(new RawDataForDataGrid(filepath));
                    }
                    break;

                default:
                    break;
            }
        }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            bool result = true;
            foreach (var character in e.Text)
            {
                if (char.IsDigit(character) && !(character == '.') && !(character == '-'))
                {
                    result = false;
                }
            }
            e.Handled = result;
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            AnalysisStatusTextBox.Clear();

            BackgroundWorker worker = new BackgroundWorker();
            worker.WorkerReportsProgress = true;
            worker.DoWork += Worker_DoWork;
            worker.ProgressChanged += Worker_ProgressChanged;
            worker.RunWorkerAsync();
        }

        void Worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            AnalysisStatusTextBox.Text = e.UserState.ToString();
        }

        private void ToggleButtons()
        {
            AnalyzeButton.IsEnabled = AllowFileDrop;
            AddDataButton.IsEnabled = AllowFileDrop;
            AddDatabaseButton.IsEnabled = AllowFileDrop;
            ClearDataButton.IsEnabled = AllowFileDrop;
            ClearDatabaseButton.IsEnabled = AllowFileDrop;
            UseBadRatiosCheckBox.IsEnabled = AllowFileDrop;
            MinValidValuesTimepointTextBox.IsEnabled = AllowFileDrop;
            MinValidValuesTotalTextBox.IsEnabled = AllowFileDrop;
            MetaMorpheusRadioButton.IsEnabled = AllowFileDrop;
            MaxQuantRadioButton.IsEnabled = AllowFileDrop;
            OutputFigures.IsEnabled = AllowFileDrop;
            RemoveBadPeptidesCheckBox.IsEnabled = AllowFileDrop;
        }

        private void Worker_DoWork(object sender, DoWorkEventArgs e)
        {
            Dispatcher.Invoke(() =>
            {
                FilesToDisplayObservableCollection.Clear();
                FilesToHideObservableCollection.Clear();
            });

            //check if we can load old results
            List<string> quantifiedPeptideInputFiles = DataFilesObservableCollection.Select(x => x.FilePath).ToList();
            List<string> turnoverResultFiles = new List<string>();
            foreach (string file in DataFilesObservableCollection.Select(x => x.FilePath))
            {
                string inputFile = file;
                string directory = Directory.GetParent(inputFile).FullName;
                string filename = Path.GetFileNameWithoutExtension(inputFile); //can be either the fast load file or the original input
                
                string fastLoadFile = Path.Combine(directory, filename + "_ApplETurnoverSavedSession.tsv");
                if (File.Exists(fastLoadFile))
                {
                    Dispatcher.Invoke(() =>
                    {
                        FilesToDisplayObservableCollection.Add(inputFile);
                    });

                    quantifiedPeptideInputFiles.Remove(inputFile); //remove so we don't analyze it again
                    quantifiedPeptideInputFiles.Remove(fastLoadFile); //remove so we don't analyze it again
                    if (!DataPreparation.LoadExistingResults(inputFile, fastLoadFile, PoolParameterDictionary, AllPeptides, AnalyzedProteins, AnalyzedProteoforms))
                    {
                        //something went wrong. Bail, reset, and run normally
                        quantifiedPeptideInputFiles = ResetFastLoadAttempt();
                        break;
                    }
                }
                else
                {
                    //something went wrong. Bail, reset, and run normally
                    quantifiedPeptideInputFiles = ResetFastLoadAttempt();
                    break;
                }
            }

            if (quantifiedPeptideInputFiles.Count != 0 && DatabasesObservableCollection.Count != 0)
            {
                AllowFileDrop = false;
                Dispatcher.Invoke(() =>
                {
                    ToggleButtons();
                    FilesToDisplayObservableCollection.Clear();
                    FilesToHideObservableCollection.Clear();
                });

                try
                {
                    (sender as BackgroundWorker).ReportProgress(0, "Starting");


                    Settings settings = GetUserSpecifiedSettings();

                    string maxStatus = quantifiedPeptideInputFiles.Count.ToString();
                    int status = 1;

                    //reading database
                    (sender as BackgroundWorker).ReportProgress(0, "Reading Database");
                    List<Protein> theoreticalProteins = DataPreparation.LoadProteins(DatabasesObservableCollection.Select(x => x.FilePath).ToList()).OrderBy(x => x.Accession).ToList();
                    Dispatcher.Invoke(() =>
                    {
                        PeptidesToDisplay.Clear();
                        AllPeptides.Clear();
                    });

                    foreach (string originalFile in quantifiedPeptideInputFiles)
                    {
                        string file = originalFile.Replace("_ApplETurnoverSavedSession", ""); //remove the extension if there was a failed load in a multi-file analysis
                        string path = Path.Combine(Directory.GetParent(file).FullName, Path.GetFileNameWithoutExtension(file));
                        ////check if the file has already been analyzed
                        //if (!File.Exists(path + "_TurnoverResults.txt"))
                        //{
                        //Load data, filter, process, parsimony
                        (sender as BackgroundWorker).ReportProgress(0, "Reading File " + status.ToString() + "/" + maxStatus + "...");
                        List<PeptideTurnoverObject> peptides = DataPreparation.ReadData(file, settings, theoreticalProteins);

                        if (peptides.Count == 0)
                        {
                            throw new Exception("No peptides were found for file: " + file + "; did you select the correct search engine?");
                        }
                            //Fit data to model, get half lives, confidence intervals
                            (sender as BackgroundWorker).ReportProgress(0, "Analyzing File " + status.ToString() + "/" + maxStatus + "...");
                        //debug
                        peptides = peptides.OrderBy(x => x.FullSequence).ToList();
                        for(int i=1; i<peptides.Count; i++)
                        {
                            if(peptides[i].FullSequence.Equals(peptides[i-1].FullSequence))
                            { }
                        }

                        PoolParameters poolParams = NonLinearRegression.RegressionAnalysis(peptides, file, settings);

                        //get protein info
                        var proteinGroups = peptides.GroupBy(x => x.Protein).ToList();
                        List<PeptideTurnoverObject> proteins = NonLinearRegression.GetProteinInfo(peptides, file, proteinGroups, "Protein");


                        //get proteoform info
                        var proteoformGroups = peptides.GroupBy(x => x.Proteoform).ToList();
                        List<PeptideTurnoverObject> proteoforms = NonLinearRegression.GetProteinInfo(peptides, file, proteoformGroups, "Proteoform");

                        AnalyzedProteins.AddRange(proteins);
                        AnalyzedProteoforms.AddRange(proteoforms);
                        PoolParameterDictionary.Add(file, poolParams);
                        PlotFit(poolParams, Path.GetFileNameWithoutExtension(file)+ " Free Amino Acids");

                        //save results to allow for quick loading in the future
                        string directory = Directory.GetParent(file).FullName;
                        string filename = Path.GetFileNameWithoutExtension(file);
                        string resultFile = Path.Combine(directory, filename + "_ApplETurnoverSavedSession.tsv");
                        DataPreparation.WriteQuickLoadFile(resultFile, poolParams, peptides, proteins, proteoforms);

                        //Add the peptides/proteins to the collection for viewing in the GUI
                        Dispatcher.Invoke(() =>
                        {
                            foreach (PeptideTurnoverObject peptide in peptides)
                            {
                                AllPeptides.Add(peptide);
                                PeptidesToDisplay.Add(peptide);
                            }
                            DisplayPeptidesDataGrid.Items.Refresh();
                            FilesToDisplayObservableCollection.Add(file);
                        });

                        status++;
                    }

                    (sender as BackgroundWorker).ReportProgress(0, "Running Statistics");
                    TTest.CompareProteinsAcrossFiles(quantifiedPeptideInputFiles, AnalyzedProteins, PoolParameterDictionary);

                    TTest.CompareProteoformsWithinFiles(quantifiedPeptideInputFiles, AnalyzedProteoforms, PoolParameterDictionary);

                    (sender as BackgroundWorker).ReportProgress(0, "Finished!");
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Task failed: " + ex.Message);
                    (sender as BackgroundWorker).ReportProgress(0, "Task failed");
                }

                AllowFileDrop = true;
                Dispatcher.Invoke(() =>
                {
                    ToggleButtons();
                });
            }
            else if (DataFilesObservableCollection.Count() != 0 && AllPeptides.Count != 0) //we had files to analyze but we were able to fast load them
            {
                Dispatcher.Invoke(() =>
                {
                    foreach (PeptideTurnoverObject peptide in AllPeptides)
                    {
                        PeptidesToDisplay.Add(peptide);
                    }
                });
                (sender as BackgroundWorker).ReportProgress(0, "Finished!");
            }
            else
            {
                MessageBox.Show("Input files are missing. The run has been stopped.");
                (sender as BackgroundWorker).ReportProgress(0, "Task failed: Input missing files.");
            }
        }

        private void ClearDataFiles_Click(object sender, RoutedEventArgs e)
        {
            DataFilesObservableCollection.Clear();
            //UpdateOutputFolderTextbox();
        }

        private void PlotFit(PoolParameters poolParams, string label, double kbi = 1000000) //default kbi is just a very high number to simulate instantaneous turnover (proxy for the available amino acid pool)
        {
            double[] timepoints = new double[1000];
            for (int i = 0; i < timepoints.Length; i++)
            {
                timepoints[i] = i / 10.0;
            }

            //half-life = ln(2)/kbt, make half life 0, kbt = infinity
            double[] rfs = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(poolParams.Kst, poolParams.Kbt, poolParams.Kao, kbi, timepoints);
            Dispatcher.Invoke(() =>
            {
                RatioComparisonPlot.plt.Layout(titleHeight: 20, xLabelHeight: 40, y2LabelWidth: 20);
                RatioComparisonPlot.plt.XLabel("Time (Days)", fontSize: 20);
                RatioComparisonPlot.plt.YLabel("Relative Fraction (Lys0/Total)", fontSize: 20);
                RatioComparisonPlot.plt.Axis(0, 100, 0, 1);
                RatioComparisonPlot.plt.Ticks(fontSize: 18);
                RatioComparisonPlot.plt.PlotScatter(timepoints, rfs, label: label, markerShape: ScottPlot.MarkerShape.none);
                if (DisplayLegendCheckBox.IsChecked.Value)
                {
                    RatioComparisonPlot.plt.Legend();
                }
                RatioComparisonPlot.Render();
            });
        }

        private Settings GetUserSpecifiedSettings()
        {
            Settings settingsToReturn = new Settings();
            Dispatcher.Invoke(() =>
            {
                settingsToReturn = new Settings(
                    int.Parse(MinValidValuesTotalTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture),
                    int.Parse(MinValidValuesTimepointTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture),
                    UseBadRatiosCheckBox.IsChecked.Value,
                    MaxQuantRadioButton.IsChecked.Value ? Settings.SearchEngine.MaxQuant : Settings.SearchEngine.MetaMorpheus,
                    RemoveBadPeptidesCheckBox.IsChecked.Value);
            });
            return settingsToReturn;
        }

        private void AddDataFile_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "MetaMorpheus Files(*.psmtsv;*.tsv;*.txt)|*.psmtsv;*.tsv;*.txt",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            }

            dataGridPeptideFiles.Items.Refresh();
        }

        private void AddDatabase_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Database(*.fasta;*.xml)|*.fasta;*.xml",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            }

            dataGridPeptideFiles.Items.Refresh();
        }

        private void ClearDatabase_Click(object sender, RoutedEventArgs e)
        {
            DatabasesObservableCollection.Clear();
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void PeptideDataGrid_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            UpdatePeptidePlots();
        }

        private void UpdatePeptidePlots()
        {
            //if nothing is selected, plot the previous selection (if any)
            if (DisplayPeptidesDataGrid.SelectedItem == null)
            {
                if (PeptidesPreviouslyPlotted.Count != 0 && PeptidesPreviouslyPlotted.All(x => PeptidesToDisplay.Contains(x)))
                {
                    PlotPeptideData(PeptidesPreviouslyPlotted);
                }
                //else do nothing
            }
            else
            {
                // plot the selected peptide/protein
                PeptideTurnoverObject row = (PeptideTurnoverObject)DisplayPeptidesDataGrid.SelectedItem;
                List<PeptideTurnoverObject> peptidesToPlot = PeptidesToDisplay.Where(x => x.FullSequence.Equals(row.FullSequence) && ((DisplayProteinInSpecificTable && x.Protein.Equals(row.Protein)) || x.Proteoform.Equals(row.Proteoform))).ToList();
                for (int i = peptidesToPlot.Count - 1; i >= 0; i--)
                {
                    if (FilesToHideObservableCollection.Contains(peptidesToPlot[i].FileName))
                    {
                        peptidesToPlot.RemoveAt(i);
                    }
                }

                PlotPeptideData(peptidesToPlot);
                PeptidesPreviouslyPlotted = peptidesToPlot;
            }
        }

        private void PlotPeptideData(List<PeptideTurnoverObject> peptidesToPlot)
        {
            RatioComparisonPlot.plt.Clear();
            RatioComparisonPlot.plt.Legend(false);
            HalfLifeComparisonPlot.plt.GetPlottables().Clear();
            if (PlotAminoAcidPoolCheckBox.IsChecked.Value)
            {
                foreach (string file in FilesToDisplayObservableCollection)
                {
                    PlotFit(PoolParameterDictionary[file], Path.GetFileNameWithoutExtension(file)+" Free Amino Acids");
                }
            }

            double minError = double.PositiveInfinity;
            double maxError = double.NegativeInfinity;
            double minHalfLife = double.PositiveInfinity;
            double maxHalfLife = double.NegativeInfinity;

            foreach (PeptideTurnoverObject peptide in peptidesToPlot)
            {
                //get the title
                int fontSize =  Math.Max(Math.Min(24, 100 / (int)Math.Round(Math.Sqrt(peptide.DisplayPeptideSequence.Length))),12);
                RatioComparisonPlot.plt.Title(peptide.DisplayPeptideSequence, fontSize: fontSize);

                string protein = peptide.DisplayProteinOrProteoform;
                string filepath = peptide.FileName;
                string filename = Path.GetFileNameWithoutExtension(filepath);
                //plot actual data
                RatioComparisonPlot.plt.PlotScatter(peptide.Timepoints, peptide.RelativeFractions, markerSize: 4, lineWidth: 0, label: filename + " Observed Ratios");

                //Plot protein info
                List<PeptideTurnoverObject> peptidesSharingProteinAndFile = PeptidesToDisplay.Where(x => x.DisplayProteinOrProteoform.Equals(protein) && x.FileName.Equals(filepath)).ToList();
                double[] errors = peptidesSharingProteinAndFile.Select(x => x.Error).ToArray();
                double[] halfLives = peptidesSharingProteinAndFile.Select(x => Math.Log(2, Math.E) / x.Kbi).ToArray();
                double[] negativeErrors = peptidesSharingProteinAndFile.Select(x => (Math.Log(2, Math.E) / x.Kbi) - (Math.Log(2, Math.E) / x.HighKbi)).ToArray();
                double[] positiveErrors = peptidesSharingProteinAndFile.Select(x => (Math.Log(2, Math.E) / x.LowKbi) - (Math.Log(2, Math.E) / x.Kbi)).ToArray();

                HalfLifeComparisonPlot.plt.Title(protein, fontSize: 24);
                HalfLifeComparisonPlot.plt.Layout(titleHeight: 20, xLabelHeight: 40, y2LabelWidth: 20);
                HalfLifeComparisonPlot.plt.YLabel("Half-life (Days)", fontSize: 20);
                HalfLifeComparisonPlot.plt.XLabel("Error (MSE)", fontSize: 20);
                HalfLifeComparisonPlot.plt.Ticks(fontSize: 18);

                double errorDiff = errors.Max() - errors.Min();
                if (errorDiff == 0)
                {
                    errorDiff = 0.01;
                }
                double halflifeDiff = halfLives.Max() - halfLives.Min();
                if (halflifeDiff == 0)
                {
                    halflifeDiff = 0.01;
                }
                var scatter = HalfLifeComparisonPlot.plt.PlotScatter(errors, halfLives, lineWidth: 0, label: filename + " peptides");
                //plot the single point of the selected peptie separately (overlay) so that we know which one it is
                var point = HalfLifeComparisonPlot.plt.PlotPoint(peptide.Error, Math.Log(2, Math.E) / peptide.Kbi);
                //plot errors
                HalfLifeComparisonPlot.plt.PlotErrorBars(errors, halfLives, null, null, positiveErrors, negativeErrors, scatter.color);
                HalfLifeComparisonPlot.plt.PlotErrorBars(new double[] { peptide.Error }, new double[] { Math.Log(2, Math.E) / peptide.Kbi },
                    null, null, new double[] { (Math.Log(2, Math.E) / peptide.LowKbi) - (Math.Log(2, Math.E) / peptide.Kbi) },
                    new double[] { (Math.Log(2, Math.E) / peptide.Kbi) - (Math.Log(2, Math.E) / peptide.HighKbi) }, color: point.color);

                minError = Math.Min(minError, errors.Min() - errorDiff * 0.2);
                maxError = Math.Max(maxError, errors.Max() + errorDiff * 0.2); 
                minHalfLife = Math.Min(minHalfLife, halfLives.Min() - negativeErrors.Max() - halflifeDiff * 0.2);
                maxHalfLife = Math.Max(maxHalfLife, halfLives.Max() + positiveErrors.Max() + halflifeDiff * 0.2);

                double ySpacingFactor = (maxHalfLife - minHalfLife) * 0.2;
                HalfLifeComparisonPlot.plt.Axis(minError, maxError, minHalfLife - ySpacingFactor, maxHalfLife + ySpacingFactor);
                HalfLifeComparisonPlot.plt.Axis();

                PeptideTurnoverObject currentProtein = DisplayProteinInSpecificTable ? 
                    AnalyzedProteins.Where(x => x.Protein.Equals(protein) && x.FileName.Equals(filepath)).FirstOrDefault() : 
                    AnalyzedProteoforms.Where(x => x.Proteoform.Equals(protein) && x.FileName.Equals(filepath)).FirstOrDefault();

                if(currentProtein==null)
                {
                    MessageBox.Show("Unable to find the protein for this peptide. There may be an issue with the loaded file.");
                    return;
                }

                //plot the fit
                if (PlotBestFitCheckBox.IsChecked.Value)
                {
                    //peptide level
                    PlotFit(PoolParameterDictionary[filepath], filename+ " Fit (" + (Math.Log(2, Math.E) / peptide.Kbi).ToString("F1") + ")", peptide.Kbi);
                    //protein level
                    double halfLife = Math.Log(2, Math.E) / currentProtein.Kbi;
                        HalfLifeComparisonPlot.plt.PlotHLine(halfLife, label: filename + " Half-life ("+halfLife.ToString("F1")+")");
                }
                //plt the confidence intervals
                if (PlotCICheckBox.IsChecked.Value)
                {
                    //peptide level
                    PlotFit(PoolParameterDictionary[filepath], filename + " Upper CI (" + (Math.Log(2, Math.E) / peptide.LowKbi).ToString("F1") + ")", peptide.LowKbi);
                    PlotFit(PoolParameterDictionary[filepath], filename + " Lower CI (" + (Math.Log(2, Math.E) / peptide.HighKbi).ToString("F1") + ")", peptide.HighKbi);
                    //protein level
                    double upperHL = Math.Log(2, Math.E) / currentProtein.LowKbi;
                    double lowerHL = Math.Log(2, Math.E) / currentProtein.HighKbi;
                        HalfLifeComparisonPlot.plt.PlotHLine(upperHL, label: filename + " Upper CI (" + upperHL.ToString("F1") + ")");
                        HalfLifeComparisonPlot.plt.PlotHLine(lowerHL, label: filename + " Lower CI (" + lowerHL.ToString("F1") + ")");
                }
            }
            if (DisplayLegendCheckBox.IsChecked.Value)
            {
                HalfLifeComparisonPlot.plt.Legend();
            }
            else
            {
                HalfLifeComparisonPlot.plt.Legend(false);
            }
            HalfLifeComparisonPlot.Render();
        }

        private void ChangeFilesToPlotSpecificData_Click(object sender, RoutedEventArgs e)
        {
            //If hiding cells
            var selectedCells = DisplayedSamplesDataGrid.SelectedCells;
            if (selectedCells != null)
            {
                foreach (var cell in selectedCells)
                {
                    string file = cell.Item.ToString();
                    FilesToDisplayObservableCollection.Remove(file);
                    FilesToHideObservableCollection.Add(file);
                }
            }

            selectedCells = HiddenSamplesDataGrid.SelectedCells;
            if (selectedCells != null)
            {
                foreach (var cell in selectedCells)
                {
                    string file = cell.Item.ToString();
                    FilesToHideObservableCollection.Remove(file);
                    FilesToDisplayObservableCollection.Add(file);
                }
            }
            UpdateDisplayedPeptides();
        }

        private void ProteinSearchTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            UpdateDisplayedPeptides();
        }

        private void UpdateDisplayedPeptides()
        {
            PeptidesToDisplay.Clear();

            //check if using search criteria
            string proteinText = ProteinSearchTextBox.Text;
            string peptideText = PeptidenSearchTextBox.Text;
            int proteinLength = proteinText.Length;
            //if filtering by both protein and peptide
            if (ProteinSearchTextBox.Text.Length >= 3 && PeptidenSearchTextBox.Text.Length >= 3)
            {
                foreach (PeptideTurnoverObject peptide in AllPeptides.Where(x => FilesToDisplayObservableCollection.Contains(x.FileName) 
                && x.Protein.Substring(0, Math.Min(proteinLength, x.Protein.Length)).Equals(proteinText) 
                && x.DisplayPeptideSequence.Contains(peptideText)))
                {
                    PeptidesToDisplay.Add(peptide);
                } 
            }
            //if filtering by just protein
            else if(ProteinSearchTextBox.Text.Length >= 3)
            {
                foreach (PeptideTurnoverObject peptide in AllPeptides.Where(x => FilesToDisplayObservableCollection.Contains(x.FileName) && x.Protein.Substring(0, Math.Min(proteinLength, x.Protein.Length)).Equals(proteinText)))
                {
                    PeptidesToDisplay.Add(peptide);
                }
            }
            //if filtering by just peptide
            else if(PeptidenSearchTextBox.Text.Length >= 3)
            {
                foreach (PeptideTurnoverObject peptide in AllPeptides.Where(x => FilesToDisplayObservableCollection.Contains(x.FileName) && x.DisplayPeptideSequence.Contains(peptideText)))
                {
                    PeptidesToDisplay.Add(peptide);
                }
            }
            //if not filtering
            else
            {
                foreach (PeptideTurnoverObject peptide in AllPeptides.Where(x => FilesToDisplayObservableCollection.Contains(x.FileName)))
                {
                    PeptidesToDisplay.Add(peptide);
                }
            }

            //remove duplicate entries (occurs if moving from proteoform to protein display)
            for(int i=0; i<PeptidesToDisplay.Count; i++)
            {
                var currentPeptide = PeptidesToDisplay[i];
                for (int j = i + 1; j < PeptidesToDisplay.Count; j++)
                {
                    var nextPeptide = PeptidesToDisplay[j];
                    if (currentPeptide.FullSequence.Equals(nextPeptide.FullSequence) && currentPeptide.FileName.Equals(nextPeptide.FileName))
                    {
                        PeptidesToDisplay.RemoveAt(j);
                        j--;
                    }
                }
            }

            DisplayedSamplesDataGrid.Items.Refresh();
            HiddenSamplesDataGrid.Items.Refresh();

            UpdatePeptidePlots();
        }

        private void PlotPrecisionScatterPlot(List<PeptideTurnoverObject> peptidesToPlot, PoolParameters poolParams)
        {
            PrecisionPlot.plt.Clear();
            PrecisionPlot.plt.GetPlottables().Clear();

            if (peptidesToPlot.Count == 0)
            {
                return;
            }

            Dictionary<double, List<(double halfLife, double relativeFraction)>> dictionaryToPlot = new Dictionary<double, List<(double halfLife, double relativeFraction)>>();

            foreach (PeptideTurnoverObject peptide in peptidesToPlot)
            {
                //grab measurements
                double halfLife = Math.Log(2, Math.E) / peptide.Kbi;
                for (int i = 0; i < peptide.Timepoints.Length; i++)
                {
                    if (dictionaryToPlot.ContainsKey(peptide.Timepoints[i]))
                    {
                        dictionaryToPlot[peptide.Timepoints[i]].Add((halfLife, peptide.RelativeFractions[i]));
                    }
                    else
                    {
                        dictionaryToPlot[peptide.Timepoints[i]] = new List<(double halfLife, double relativeFraction)> { (halfLife, peptide.RelativeFractions[i]) };
                    }
                }
            }

            //plot all peptide data
            double[] timepoints = dictionaryToPlot.Keys.OrderBy(x => x).ToArray();
            foreach (double timepoint in timepoints)
            {
                var value = dictionaryToPlot[timepoint];
                PrecisionPlot.plt.PlotScatter(value.Select(x => x.halfLife).ToArray(), value.Select(x => x.relativeFraction).ToArray(), lineWidth: 0, markerSize: 3, label: timepoint.ToString(), markerShape: ScottPlot.MarkerShape.openCircle);
            }

            //plt fits for each timepoint on top of the peptide data
            double[] halflives = new double[2499];
            for (int i = 0; i < halflives.Length; i++)
            {
                halflives[i] = i / 5.0 + 0.2;
            }

            List<double>[] rfs = new List<double>[timepoints.Length];
            for (int i = 0; i < timepoints.Length; i++)
            {
                rfs[i] = new List<double>();
            }

            foreach (double halflife in halflives)
            {
                if (halflife == 0)
                {
                    continue;
                }

                //halflife = ln(2)/kbi
                //kbi = ln(2)/halflife


                double[] rfsForThisHalfLife = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(
                    poolParams.Kst, poolParams.Kbt, poolParams.Kao, Math.Log(2, Math.E) / halflife, timepoints);

                for (int i = 0; i < timepoints.Length; i++)
                {
                    rfs[i].Add(rfsForThisHalfLife[i]);
                }
            }
            PrecisionPlot.plt.Axis(0, 50, 0, 1);
            for (int i = 0; i < timepoints.Length; i++)
            {
                PrecisionPlot.plt.PlotScatter(halflives, rfs[i].ToArray(), Color.Black, markerSize: 0);
            }

            if (DisplayLegendCheckBox.IsChecked.Value)
            {
                PrecisionPlot.plt.Legend(location: ScottPlot.legendLocation.upperRight);
            }
            else
            {
                PrecisionPlot.plt.Legend(false);
            }

            PrecisionPlot.plt.YLabel("Lys0/Total");
            PrecisionPlot.plt.XLabel("Half-life (Days)");
            PrecisionPlot.plt.Axis(0, 50, 0, 1);
            PrecisionPlot.Render(); PrecisionPlot.Render();
        }

        private void PlotHalfLifeHistogram(List<PeptideTurnoverObject> peptidesToPlot)
        {
            HalfLifeHistogramPlot.plt.Clear();
            HalfLifeHistogramPlot.plt.GetPlottables().Clear();
            double[] halflives = peptidesToPlot.Select(x => Math.Log(2, Math.E) / x.Kbi).ToArray();
            Histogram histogram = new Histogram(halflives, min: 0, max: 100);

            double barWidth = histogram.binSize * 1.2; // slightly over-side to reduce anti-alias rendering artifacts

            HalfLifeHistogramPlot.plt.Axis(0, 100, 0, histogram.counts.Max() * 1.1);
            HalfLifeHistogramPlot.plt.PlotBar(histogram.bins, histogram.counts, barWidth: barWidth, outlineWidth: 0);
            HalfLifeHistogramPlot.plt.YLabel("Frequency (# peptides)");
            HalfLifeHistogramPlot.plt.XLabel("Half-life (Days)");
            HalfLifeHistogramPlot.plt.Axis(0, 100, 0, histogram.counts.Max() * 1.1);
            HalfLifeComparisonPlot.Render(); HalfLifeComparisonPlot.Render();
        }

        private void ExportFiguresButton_Click(object sender, RoutedEventArgs e)
        {
            //create a new directory to export the figures
            try
            {
                var pathOfFirstSpectraFile = Path.GetDirectoryName(DataFilesObservableCollection.First().FilePath);

                var exportFolder = Path.Combine(pathOfFirstSpectraFile, "ExportedFigures");
                if (!Directory.Exists(exportFolder))
                {
                    Directory.CreateDirectory(exportFolder);
                }
                var saveFolder = Path.Combine(exportFolder, DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture)); //make a unique folder for each save instance
                Directory.CreateDirectory(saveFolder);

                //save the plots within the directory
                PrecisionPlot.plt.SaveFig(Path.Combine(saveFolder, "GlobalPrecisionScatterPlot.tif"));
                HalfLifeHistogramPlot.plt.SaveFig(Path.Combine(saveFolder, "GlobalHalfLifeHistogram.tif"));
                RatioComparisonPlot.plt.SaveFig(Path.Combine(saveFolder, "SinglePeptideScatterPlot.tif"));
                HalfLifeComparisonPlot.plt.SaveFig(Path.Combine(saveFolder, "SingleProteinScatterPlot.tif"));
            }
            catch
            {
                MessageBox.Show("No results found. Please run an experiment.");
            }
        }

        private void AnalyzedFilesGrid_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            UpdateGlobalVisualization();
        }

        private void GlobalVisualization_Click(object sender, RoutedEventArgs e)
        {
            UpdateGlobalVisualization();
        }


        private void UpdateGlobalVisualization(PoolParameters customParams = null)
        {
            if (DisplayAnalyzedFilesDataGrid.SelectedItem == null)
            {
                return;
            }
            //plot the precision scatter plot and the half-life histogram for this file
            string dataFile = ((RawDataForDataGrid)DisplayAnalyzedFilesDataGrid.SelectedItem).FilePath;
            if (peptideRadioButton.IsChecked.Value)
            {
                List<PeptideTurnoverObject> peptidesForThisFile = AllPeptides.Where(x => x.FileName.Equals(dataFile)).ToList();
                if (customParams != null || PoolParameterDictionary.ContainsKey(dataFile)) //if the file has been analyzed
                {
                    bool firstTime = customParams == null;
                    customParams = customParams ?? PoolParameterDictionary[dataFile];
                    ChangeParamTextBox = false;
                    if (firstTime)
                    {
                        KstTB.Text = customParams.Kst.ToString();
                        KbtTB.Text = customParams.Kbt.ToString();
                        KaoTB.Text = customParams.Kao.ToString();
                    }
                    MseTB.Text = peptidesForThisFile.Average(x => x.Error).ToString();
                    ChangeParamTextBox = true;
                    PlotPrecisionScatterPlot(peptidesForThisFile, customParams);
                    PlotHalfLifeHistogram(peptidesForThisFile);
                }
            }
            else //graph proteins
            {
                List<PeptideTurnoverObject> proteinsForThisFile = AnalyzedProteins.Where(x => x.FileName.Equals(dataFile)).ToList();
                if (customParams != null || PoolParameterDictionary.ContainsKey(dataFile)) //if the file has been analyzed
                {
                    customParams = customParams ?? PoolParameterDictionary[dataFile];
                    PlotPrecisionScatterPlot(proteinsForThisFile, customParams);
                    PlotHalfLifeHistogram(proteinsForThisFile);
                }
            }
        }

        private void ParamTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            if (ChangeParamTextBox)
            {
                try
                {
                    UpdateGlobalVisualization(new PoolParameters(Convert.ToDouble(KstTB.Text), Convert.ToDouble(KbtTB.Text), Convert.ToDouble(KaoTB.Text)));
                }
                catch
                {
                    //do nothing, the user is still typing
                }
            }
        }

        private void ParamApply_Click(object sender, RoutedEventArgs e)
        {
            if (DisplayAnalyzedFilesDataGrid.SelectedItem == null)
            {
                return;
            }
            assignParamsButton.IsEnabled = false;
            string dataFile = ((RawDataForDataGrid)DisplayAnalyzedFilesDataGrid.SelectedItem).FilePath;
            PoolParameters customParams = new PoolParameters(Convert.ToDouble(KstTB.Text), Convert.ToDouble(KbtTB.Text), Convert.ToDouble(KaoTB.Text));
            PoolParameterDictionary[dataFile] = customParams;
            List<PeptideTurnoverObject> peptidesForThisFile = AllPeptides.Where(x => x.FileName.Equals(dataFile)).ToList();

            NonLinearRegression.UpdateKbi(customParams.Kst, customParams.Kbt, customParams.Kao, peptidesForThisFile, 0.001);
            NonLinearRegression.UpdateKbi(customParams.Kst, customParams.Kbt, customParams.Kao, peptidesForThisFile, 0.0001);
            NonLinearRegression.UpdateKbi(customParams.Kst, customParams.Kbt, customParams.Kao, peptidesForThisFile, 0.00001);
            MseTB.Text = peptidesForThisFile.Sum(x => x.Error).ToString();
            UpdateGlobalVisualization();
            assignParamsButton.IsEnabled = true;
            //CreateMapForLocalMinimaSearch();
        }

        private void SpecificVisualization_Click(object sender, RoutedEventArgs e)
        {
            UpdateSpecificVisualization();
        }

        private void UpdateSpecificVisualization()
        {
            //change display of protein vs proteoform naming
            DisplayProteinInSpecificTable = proteinSpecificRadioButton.IsChecked.Value;
            DisplayFullPeptideSequence = showFullSequenceCheckbox.IsChecked.Value;

            //add/remove redundant peptide entries from table
            PeptidesToDisplay.Clear();
            if (DisplayProteinInSpecificTable)
            {
                string previousFullSequence = "";
                foreach (PeptideTurnoverObject peptide in AllPeptides.OrderBy(x=>x.FullSequence))
                {
              //      if (!previousFullSequence.Equals(peptide.FullSequence))
                    {
                        previousFullSequence = peptide.FullSequence;
                        peptide.UpdateDisplayProteinOrProteoform(DisplayProteinInSpecificTable, DisplayFullPeptideSequence);
                        PeptidesToDisplay.Add(peptide);
                    }
                }
            }
            else
            {
                foreach (PeptideTurnoverObject peptide in AllPeptides)
                {
                    peptide.UpdateDisplayProteinOrProteoform(DisplayProteinInSpecificTable, DisplayFullPeptideSequence);
                    PeptidesToDisplay.Add(peptide);
                }
            }
            UpdateDisplayedPeptides();

            //refresh
            DisplayedSamplesDataGrid.Items.Refresh();
            HiddenSamplesDataGrid.Items.Refresh();

       //     UpdatePeptidePlots();
        }

        private void AddRemoveLegends_Click(object sender, RoutedEventArgs e)
        {
            ChangeFilesToPlotSpecificData_Click(sender, e);
            UpdateGlobalVisualization();
        }

        private List<string> ResetFastLoadAttempt()
        {            
            PoolParameterDictionary.Clear();
            Dispatcher.Invoke(() =>
            {
                FilesToDisplayObservableCollection.Clear();
            });
            AllPeptides.Clear();
            AnalyzedProteins.Clear();
            AnalyzedProteoforms.Clear();
            return DataFilesObservableCollection.Select(x => x.FilePath).ToList();
        }
    }
}