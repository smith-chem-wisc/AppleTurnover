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
        private List<PeptideTurnoverObject> PeptidesPreviouslyPlotted = new List<PeptideTurnoverObject>();
        private bool ChangeParamTextBox;

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
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (AllowFileDrop)
            {
                string[] files = ((string[])e.Data.GetData(DataFormats.FileDrop)).OrderBy(p => p).ToArray();

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
                string filename = Path.GetFileNameWithoutExtension(inputFile);
                //check the user didn't grab the wrong file
                if (filename.Contains("_ApplETurnoverSavedSession"))
                {
                    filename = filename.Replace("_ApplETurnoverSavedSession", "");
                    string updatedPath = inputFile.Replace("_ApplETurnoverSavedSession", "");
                    //if the user grabbed the saved session and the original input file, remove the session
                    var badFile = DataFilesObservableCollection.Where(x => x.FilePath.Equals(inputFile)).FirstOrDefault();

                    if (DataFilesObservableCollection.Select(x => x.FilePath).Contains(updatedPath))
                    {
                        DataFilesObservableCollection.Remove(badFile);
                    }
                    else //treat this as the original input file
                    {
                        badFile.RemoveSessionTag();
                    }
                    inputFile = updatedPath;
                    Dispatcher.Invoke(() =>
                    {
                        dataGridPeptideFiles.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridPeptideFiles.Items.Refresh();
                    });
                }
                string resultFile = Path.Combine(directory, filename + "_ApplETurnoverSavedSession.tsv");
                if (File.Exists(resultFile))
                {
                    Dispatcher.Invoke(() =>
                    {
                        FilesToDisplayObservableCollection.Add(inputFile);
                    });

                    quantifiedPeptideInputFiles.Remove(inputFile); //remove so we don't analyze it again
                    if (!DataPreparation.LoadExistingResults(inputFile, resultFile, PoolParameterDictionary, AllPeptides, AnalyzedProteins))
                    {
                        //something went wrong. Bail, reset, and run normally
                        quantifiedPeptideInputFiles = DataFilesObservableCollection.Select(x => x.FilePath).ToList();
                        PoolParameterDictionary.Clear();
                        FilesToDisplayObservableCollection.Clear();
                        AllPeptides.Clear();
                        AnalyzedProteins.Clear();
                        break;
                    }
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

                // try
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

                    foreach (string file in quantifiedPeptideInputFiles)
                    {
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
                        PoolParameters poolParams = NonLinearRegression.RegressionAnalysis(peptides, file, settings);
                        //get protein info
                        List<PeptideTurnoverObject> proteins = NonLinearRegression.GetProteinInfo(peptides, file);
                        AnalyzedProteins.AddRange(proteins);
                        PoolParameterDictionary.Add(file, poolParams);
                        PlotFit(poolParams, "Free Amino Acids");

                        //save results to allow for quick loading in the future
                        string directory = Directory.GetParent(file).FullName;
                        string filename = Path.GetFileNameWithoutExtension(file);
                        string resultFile = Path.Combine(directory, filename + "_ApplETurnoverSavedSession.tsv");
                        DataPreparation.WriteResults(resultFile, poolParams, peptides, proteins);

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

                    TTest.CompareProteoformsWithinFiles(quantifiedPeptideInputFiles, AnalyzedProteins, PoolParameterDictionary);

                    (sender as BackgroundWorker).ReportProgress(0, "Finished!");
                }
                //catch (Exception ex)
                //{
                //    MessageBox.Show("Task failed: " + ex.Message);
                //    (sender as BackgroundWorker).ReportProgress(0, "Task failed");
                //}

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
            double[] timepoints = new double[200];
            for (int i = 0; i < timepoints.Length; i++)
            {
                timepoints[i] = i / 2.0;
            }

            //half-life = ln(2)/kbt, make half life 0, kbt = infinity
            double[] rfs = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(poolParams.Kst, poolParams.Kbt, poolParams.Kao, kbi, timepoints);
            Dispatcher.Invoke(() =>
            {
                RatioComparisonPlot.plt.Layout(titleHeight: 20, xLabelHeight: 40, y2LabelWidth: 20);
                RatioComparisonPlot.plt.XLabel("Time (Days)", fontSize: 24);
                RatioComparisonPlot.plt.YLabel("Relative Fraction (New/Total)", fontSize: 24);
                RatioComparisonPlot.plt.Axis(0, 100, 0, 1);
                RatioComparisonPlot.plt.Ticks(fontSize: 18);
                RatioComparisonPlot.plt.PlotScatter(timepoints, rfs, label: label, markerShape: ScottPlot.MarkerShape.none);
                //      RatioComparisonPlot.plt.Legend(fontSize: 8);
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
                List<PeptideTurnoverObject> peptidesToPlot = PeptidesToDisplay.Where(x => x.FullSequence.Equals(row.FullSequence) && x.Protein.Equals(row.Protein)).ToList();
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
                    PlotFit(PoolParameterDictionary[file], "Free Amino Acids");
                }
            }

            double[] timepoints = new double[200];
            for (int i = 0; i < timepoints.Length; i++)
            {
                timepoints[i] = i / 2.0;
            }

            foreach (PeptideTurnoverObject peptide in peptidesToPlot)
            {
                //get the title
                RatioComparisonPlot.plt.Title(peptide.FullSequence, fontSize: 24);
                //plot actual data
                RatioComparisonPlot.plt.PlotScatter(peptide.Timepoints, peptide.RelativeFractions, markerSize: 4, lineWidth: 0, label: "Observed Ratios");

                //Plot protein info
                string protein = peptide.Protein;
                string file = peptide.FileName;
                List<PeptideTurnoverObject> peptidesSharingProteinAndFile = PeptidesToDisplay.Where(x => x.Protein.Equals(protein) && x.FileName.Equals(file)).ToList();
                double[] errors = peptidesSharingProteinAndFile.Select(x => x.Error).ToArray();
                double[] halfLives = peptidesSharingProteinAndFile.Select(x => Math.Log(2, Math.E) / x.Kbi).ToArray();
                double[] negativeErrors = peptidesSharingProteinAndFile.Select(x => (Math.Log(2, Math.E) / x.Kbi) - (Math.Log(2, Math.E) / x.HighKbi)).ToArray();
                double[] positiveErrors = peptidesSharingProteinAndFile.Select(x => (Math.Log(2, Math.E) / x.LowKbi) - (Math.Log(2, Math.E) / x.Kbi)).ToArray();
                HalfLifeComparisonPlot.plt.Title(protein, fontSize: 24);
                HalfLifeComparisonPlot.plt.Layout(titleHeight: 20, xLabelHeight: 40, y2LabelWidth: 20);
                HalfLifeComparisonPlot.plt.YLabel("Half-life (Days)", fontSize: 24);
                HalfLifeComparisonPlot.plt.XLabel("Error (MSE)", fontSize: 24);
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
                HalfLifeComparisonPlot.plt.Axis(errors.Min() - errorDiff * 0.2, errors.Max() + errorDiff * 0.2, halfLives.Min() - negativeErrors.Max() - halflifeDiff * 0.2, halfLives.Max() + positiveErrors.Max() + halflifeDiff * 0.2);
                var scatter = HalfLifeComparisonPlot.plt.PlotScatter(errors, halfLives, lineWidth: 0);
                //plot the single point of the selected peptie separately (overlay) so that we know which one it is
                var point = HalfLifeComparisonPlot.plt.PlotPoint(peptide.Error, Math.Log(2, Math.E) / peptide.Kbi);
                //plot errors
                HalfLifeComparisonPlot.plt.PlotErrorBars(errors, halfLives, null, null, positiveErrors, negativeErrors, scatter.color);
                HalfLifeComparisonPlot.plt.PlotErrorBars(new double[] { peptide.Error }, new double[] { Math.Log(2, Math.E) / peptide.Kbi },
                    null, null, new double[] { (Math.Log(2, Math.E) / peptide.LowKbi) - (Math.Log(2, Math.E) / peptide.Kbi) },
                    new double[] { (Math.Log(2, Math.E) / peptide.Kbi) - (Math.Log(2, Math.E) / peptide.HighKbi) }, color: point.color);

                HalfLifeComparisonPlot.plt.Axis(errors.Min() - errorDiff * 0.2, errors.Max() + errorDiff * 0.2, halfLives.Min() - negativeErrors.Max() - halflifeDiff * 0.2, halfLives.Max() + positiveErrors.Max() + halflifeDiff * 0.2);
                HalfLifeComparisonPlot.plt.Axis();

                PeptideTurnoverObject currentProtein = AnalyzedProteins.Where(x => x.Protein.Equals(protein) && x.FileName.Equals(peptide.FileName)).FirstOrDefault();

                //plot the fit
                if (PlotBestFitCheckBox.IsChecked.Value)
                {
                    //peptide level
                    PlotFit(PoolParameterDictionary[peptide.FileName], "Fit (MSE: " + peptide.Error.ToString("0.0E0") + ")", peptide.Kbi);
                    //protein level
                    HalfLifeComparisonPlot.plt.PlotHLine(Math.Log(2, Math.E) / currentProtein.Kbi);
                }
                //plt the confidence intervals
                if (PlotCICheckBox.IsChecked.Value)
                {
                    //peptide level
                    PlotFit(PoolParameterDictionary[peptide.FileName], ("Upper CI"), peptide.LowKbi);
                    PlotFit(PoolParameterDictionary[peptide.FileName], ("Lower CI"), peptide.HighKbi);
                    //protein level
                    HalfLifeComparisonPlot.plt.PlotHLine(Math.Log(2, Math.E) / currentProtein.LowKbi);
                    HalfLifeComparisonPlot.plt.PlotHLine(Math.Log(2, Math.E) / currentProtein.HighKbi);
                }

                HalfLifeComparisonPlot.Render();
            }
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

            //if using search criteria
            if (ProteinSearchTextBox.Text.Length >= 3)
            {
                string text = ProteinSearchTextBox.Text;
                int length = text.Length;
                foreach (PeptideTurnoverObject peptide in AllPeptides.Where(x => FilesToDisplayObservableCollection.Contains(x.FileName) && x.Protein.Substring(0, Math.Min(length, x.Protein.Length)).Equals(text)))
                {
                    PeptidesToDisplay.Add(peptide);
                }
            }
            else
            {
                foreach (PeptideTurnoverObject peptide in AllPeptides.Where(x => FilesToDisplayObservableCollection.Contains(x.FileName)))
                {
                    PeptidesToDisplay.Add(peptide);
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
                PrecisionPlot.plt.PlotScatter(value.Select(x => x.halfLife).ToArray(), value.Select(x => x.relativeFraction).ToArray(), lineWidth: 0, markerSize: 1.5, label: timepoint.ToString());
            }

            //plt fits for each timepoint on top of the peptide data
            double[] halflives = new double[499];
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
            //RatioComparisonPlot.plt.Layout(xLabelHeight: 40, y2LabelWidth: 20);
            //PrecisionPlot.plt.Legend();
            PrecisionPlot.plt.YLabel("New/Total");
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
                    customParams = customParams ?? PoolParameterDictionary[dataFile];
                    ChangeParamTextBox = false;
                    KstTB.Text = customParams.Kst.ToString();
                    KbtTB.Text = customParams.Kbt.ToString();
                    KaoTB.Text = customParams.Kao.ToString();
                    MseTB.Text = peptidesForThisFile.Sum(x => x.Error).ToString();
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

        private void CreateMapForLocalMinimaSearch()
        {

            //get original params
            var ogParams = PoolParameterDictionary.First().Value;
            double kst = ogParams.Kst;
            double kbt = ogParams.Kbt;
            double kao = ogParams.Kao;

            string dataFile = ((RawDataForDataGrid)DisplayAnalyzedFilesDataGrid.SelectedItem).FilePath;
            PoolParameters customParams = new PoolParameters(Convert.ToDouble(KstTB.Text), Convert.ToDouble(KbtTB.Text), Convert.ToDouble(KaoTB.Text));
            PoolParameterDictionary[dataFile] = customParams;
            List<PeptideTurnoverObject> peptidesForThisFile = AllPeptides.Where(x => x.FileName.Equals(dataFile)).ToList();

            List<string> results = new List<string>();

            double[] ratiosForIteration = new double[] { 0.1, 0.2, 0.33, 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 2, 4 };
            for (int i = 0; i < ratiosForIteration.Length; i++)
            {
                double kstCurrent = kst * ratiosForIteration[i];
                for (int j = 0; j < ratiosForIteration.Length; j++)
                {
                    double kbtCurrent = kbt * ratiosForIteration[j];
                    for (int k = 0; k < ratiosForIteration.Length; k++)
                    {
                        double kaoCurrent = kao * ratiosForIteration[k];

                        NonLinearRegression.UpdateKbi(kstCurrent, kbtCurrent, kaoCurrent, peptidesForThisFile, 0.001);
                        NonLinearRegression.UpdateKbi(kstCurrent, kbtCurrent, kaoCurrent, peptidesForThisFile, 0.0001);
                        NonLinearRegression.UpdateKbi(kstCurrent, kbtCurrent, kaoCurrent, peptidesForThisFile, 0.00001);

                        results.Add(kstCurrent.ToString() + '\t' + kbtCurrent.ToString() + '\t' + kaoCurrent.ToString() + '\t' + peptidesForThisFile.Sum(x => x.Error).ToString());
                    }
                }
            }
            File.WriteAllLines(@"E:\Chemistry\MusSILAC\TissuePaper\Grid.tsv", results);

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
            string dataFile = ((RawDataForDataGrid)DisplayAnalyzedFilesDataGrid.SelectedItem).FilePath;
            PoolParameters customParams = new PoolParameters(Convert.ToDouble(KstTB.Text), Convert.ToDouble(KbtTB.Text), Convert.ToDouble(KaoTB.Text));
            PoolParameterDictionary[dataFile] = customParams;
            List<PeptideTurnoverObject> peptidesForThisFile = AllPeptides.Where(x => x.FileName.Equals(dataFile)).ToList();

            NonLinearRegression.UpdateKbi(customParams.Kst, customParams.Kbt, customParams.Kao, peptidesForThisFile, 0.001);
            NonLinearRegression.UpdateKbi(customParams.Kst, customParams.Kbt, customParams.Kao, peptidesForThisFile, 0.0001);
            NonLinearRegression.UpdateKbi(customParams.Kst, customParams.Kbt, customParams.Kao, peptidesForThisFile, 0.00001);
            MseTB.Text = peptidesForThisFile.Sum(x => x.Error).ToString();
            UpdateGlobalVisualization();
            //CreateMapForLocalMinimaSearch();
        }
    }
}