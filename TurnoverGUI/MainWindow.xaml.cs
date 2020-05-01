using Proteomics;
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
        private Dictionary<string, PoolParameters> PoolParameterDictionary = new Dictionary<string,PoolParameters>();

        public MainWindow()
        {
            InitializeComponent();
            dataGridPeptideFiles.DataContext = DataFilesObservableCollection;
            dataGridDatabaseFiles.DataContext = DatabasesObservableCollection;
            DisplayedSamplesDataGrid.DataContext = FilesToDisplayObservableCollection;
            HiddenSamplesDataGrid.DataContext = FilesToHideObservableCollection;

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
                    if (!DataFilesObservableCollection.Any(x=> x.FilePath.Equals(filepath)))
                    {
                        DataFilesObservableCollection.Add(new RawDataForDataGrid(filepath));
                    }

                    //if (string.IsNullOrWhiteSpace(OutputFolderTextbox.Text))
                    //{
                    //    OutputFolderTextbox.Text = Path.Combine(Path.GetDirectoryName(filepath), @"$DATETIME");
                    //}

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
        }

        private void Worker_DoWork(object sender, DoWorkEventArgs e)
        {
            if (DataFilesObservableCollection.Count != 0 || DatabasesObservableCollection.Count!=0)
            {
                AllowFileDrop = false;
                Dispatcher.Invoke(() =>
                {
                    ToggleButtons();
                    FilesToDisplayObservableCollection.Clear();
                    FilesToHideObservableCollection.Clear();
                });

                //try
                {
                    (sender as BackgroundWorker).ReportProgress(0, "Starting");
                    //var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

                    //Dispatcher.Invoke(() =>
                    //{
                    //    string outputFolder = OutputFolderTextbox.Text.Replace("$DATETIME", startTimeForAllFilenames);
                    //    OutputFolderTextbox.Text = outputFolder;

                    //    if (!Directory.Exists(outputFolder))
                    //    {
                    //        Directory.CreateDirectory(outputFolder);
                    //    }
                    //});

                    Settings settings = GetUserSpecifiedSettings();

                    //TODO: Can't use protein.tsv, has messed up ratios?
                    //Take median ratio (most intense of the two if tie breaker) (What does MaxQuant do?)
                    //Flag peptides that are modified and not
                    //Can't just use base sequence, need to use modification residue
                    //have to read in the database to do this
                    //compare across all proteins for possible shared proteins
                    //interesting to keep track of this info
                    //do peptides across samples (of the same time point) correlate?
                    //if reading in multiple, how to determine significant differences?

                    string maxStatus = DataFilesObservableCollection.Count.ToString();
                    int status = 1;

                    //reading database
                    (sender as BackgroundWorker).ReportProgress(0, "Reading Database");
                    List<Protein> theoreticalProteins = DataPreparation.LoadProteins(DatabasesObservableCollection.Select(x => x.FilePath).ToList()).OrderBy(x => x.Accession).ToList();
                    Dispatcher.Invoke(() =>
                    {
                        PeptidesToDisplay.Clear();
                        AllPeptides.Clear();
                    });

                    foreach (string file in DataFilesObservableCollection.Select(x => x.FilePath))
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
                            throw new Exception("No peptides were found for file: " + file);
                        }
                            //Fit data to model, get half lives, confidence intervals
                            (sender as BackgroundWorker).ReportProgress(0, "Analyzing File " + status.ToString() + "/" + maxStatus + "...");
                        PoolParameters poolParams = NonLinearRegression.RegressionAnalysis(peptides, file);
                        PoolParameterDictionary.Add(file, poolParams);
                        PlotFit(poolParams, "Free Amino Acids");

                        //Add the peptides/proteins to the collection for viewing in the GUI
                        Dispatcher.Invoke(() =>
                        {
                            foreach (PeptideTurnoverObject peptide in peptides)
                            {
                                peptide.SetErrorString();
                                AllPeptides.Add(peptide);
                                PeptidesToDisplay.Add(peptide);
                            }
                            DisplayPeptidesDataGrid.Items.Refresh();
                            FilesToDisplayObservableCollection.Add(file);
                        });

                        status++;
                    }

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

        //private void UpdateOutputFolderTextbox()
        //{
        //    if (DataFilesObservableCollection.Any())
        //    {
        //        // if current output folder is blank and there is a spectra file, use the spectra file's path as the output path
        //        if (string.IsNullOrWhiteSpace(OutputFolderTextbox.Text))
        //        {
        //            var pathOfFirstSpectraFile = Path.GetDirectoryName(DataFilesObservableCollection.First().FilePath);
        //            OutputFolderTextbox.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
        //        }
        //        // else do nothing (do not override if there is a path already there; might clear user-defined path)
        //    }
        //    else
        //    {
        //        // no spectra files; clear the output folder from the GUI
        //        OutputFolderTextbox.Clear();
        //    }
        //}

        private void PlotFit(PoolParameters poolParams, string label, double kbi=1000000) //default kbi is just a very high number to simulate instantaneous turnover (proxy for the available amino acid pool)
        {
            double[] timepoints = new double[200];
            for(int i=0; i<timepoints.Length; i++)
            {
                timepoints[i] = i / 2.0; 
            }

            //half-life = ln(2)/kbt, make half life 0, kbt = infinity
            double[] rfs = NonLinearRegression.PredictRelativeFractionUsingThreeCompartmentModel(poolParams.Kst, poolParams.Kbt, poolParams.Koa, kbi, timepoints);
            Dispatcher.Invoke(() =>
            {
                RatioComparisonPlot.plt.Layout(titleHeight: 20, xLabelHeight: 40, y2LabelWidth: 20);
                RatioComparisonPlot.plt.XLabel("Time (Days)", fontSize: 24);
                RatioComparisonPlot.plt.YLabel("Relative Fraction (New/Total)", fontSize:24);
                RatioComparisonPlot.plt.Axis(0, 100, 0, 1);
                RatioComparisonPlot.plt.Ticks(fontSize: 18);
                RatioComparisonPlot.plt.PlotScatter(timepoints, rfs, label:label, markerShape:ScottPlot.MarkerShape.none);
                RatioComparisonPlot.plt.Legend(fontSize:8);
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
                    MaxQuantRadioButton.IsChecked.Value ? Settings.SearchEngine.MaxQuant : Settings.SearchEngine.MetaMorpheus);
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
            UpdatePlot();
        }

        private void UpdatePlot()
        {
            if (DisplayPeptidesDataGrid.SelectedItem == null)
            {
                return;
            }

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
        }

        private void PlotPeptideData(List<PeptideTurnoverObject> peptidesToPlot)
        {
            RatioComparisonPlot.plt.Clear();
            HalfLifeComparisonPlot.plt.GetPlottables().Clear();
            if (PlotAminoAcidPoolCheckBox.IsChecked.Value)
            {
                foreach(string file in FilesToDisplayObservableCollection)
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
                RatioComparisonPlot.plt.Title(peptide.FullSequence, fontSize:24);
                //plot actual data
                RatioComparisonPlot.plt.PlotScatter(peptide.Timepoints, peptide.RelativeFractions, lineWidth: 0, label: "Observed Ratios");
                //plot the fit
                PlotFit(PoolParameterDictionary[peptide.FileName], "Fit (MSE: "+peptide.Error.ToString("0.0E0")+")", peptide.Kbi);
                //plt the confidence intervals
                if (PlotCICheckBox.IsChecked.Value)
                {
                    PlotFit(PoolParameterDictionary[peptide.FileName], ("Upper CI"), peptide.LowKbi);
                    PlotFit(PoolParameterDictionary[peptide.FileName], ("Lower CI"), peptide.HighKbi);
                }

                //Plot protein info
                string protein = peptide.Protein;
                string file = peptide.FileName;
                List<PeptideTurnoverObject> peptidesSharingProteinAndFile = PeptidesToDisplay.Where(x => x.Protein.Equals(protein) && x.FileName.Equals(file)).ToList();
                double[] errors = peptidesSharingProteinAndFile.Select(x => x.Error).ToArray();
                double[] halfLives = peptidesSharingProteinAndFile.Select(x => Math.Log(2, Math.E) / x.Kbi).ToArray();
                double[] negativeErrors = peptidesSharingProteinAndFile.Select(x => (Math.Log(2, Math.E) / x.Kbi) - (Math.Log(2, Math.E) / x.HighKbi)).ToArray();
                double[] positiveErrors = peptidesSharingProteinAndFile.Select(x => (Math.Log(2, Math.E) / x.LowKbi) - (Math.Log(2, Math.E) / x.Kbi)).ToArray();
                HalfLifeComparisonPlot.plt.Title(protein, fontSize:24);
                HalfLifeComparisonPlot.plt.Layout(titleHeight: 20, xLabelHeight: 40, y2LabelWidth: 20);
                HalfLifeComparisonPlot.plt.YLabel("Half-life (Days)", fontSize: 24);
                HalfLifeComparisonPlot.plt.XLabel("Error (MSE)", fontSize: 24);
                HalfLifeComparisonPlot.plt.Ticks(fontSize: 18);
                var scatter = HalfLifeComparisonPlot.plt.PlotScatter(errors, halfLives, lineWidth:0);
                //plot the single point of the selected peptie separately (overlay) so that we know which one it is
                var point = HalfLifeComparisonPlot.plt.PlotPoint(peptide.Error, Math.Log(2, Math.E) / peptide.Kbi);
                //plot errors
                HalfLifeComparisonPlot.plt.PlotErrorBars(errors, halfLives, null, null, positiveErrors, negativeErrors, scatter.color);          
                HalfLifeComparisonPlot.plt.PlotErrorBars(new double[] { peptide.Error }, new double[] { Math.Log(2, Math.E) / peptide.Kbi },
                    null, null, new double[] { (Math.Log(2, Math.E) / peptide.LowKbi) - (Math.Log(2, Math.E) / peptide.Kbi) },
                    new double[] { (Math.Log(2, Math.E) / peptide.Kbi) - (Math.Log(2, Math.E) / peptide.HighKbi) }, color:point.color);
                HalfLifeComparisonPlot.plt.Axis(errors.Min()-0.1, errors.Max()+0.1,halfLives.Min()-negativeErrors.Max()-0.1, halfLives.Max()+positiveErrors.Max()+0.1);
                HalfLifeComparisonPlot.Render();
            }
        }

        private void ChangeFilesToPlot_Click(object sender, RoutedEventArgs e)
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
                foreach(PeptideTurnoverObject peptide in AllPeptides.Where(x => FilesToDisplayObservableCollection.Contains(x.FileName) && x.Protein.Substring(0, Math.Min(length, x.Protein.Length)).Equals(text)))
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

            UpdatePlot();
        }

        //private void OpenOutputFolder_Click(object sender, RoutedEventArgs e)
        //{
        //    string outputFolder = OutputFolderTextbox.Text;
        //    if (outputFolder.Contains("$DATETIME"))
        //    {
        //        // the exact file path isn't known, so just open the parent directory
        //        outputFolder = Directory.GetParent(outputFolder).FullName;
        //    }

        //    if (!Directory.Exists(outputFolder) && !string.IsNullOrEmpty(outputFolder))
        //    {
        //        // create the directory if it doesn't exist yet
        //            Directory.CreateDirectory(outputFolder);
        //    }

        //    if (Directory.Exists(outputFolder))
        //    {
        //        // open the directory
        //        System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo()
        //        {
        //            FileName = outputFolder,
        //            UseShellExecute = true,
        //            Verb = "open"
        //        });
        //    }
        //}
    }
}
