/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.panels;

import ca.tonita.jawbreaker.datasets.TOVFamilyDataset;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.models.JawBreakerModel;
import ca.tonita.physics.gr.hydro.TOVData;
import ca.tonita.jawbreaker.models.TOVFamily;
import ca.tonita.jawbreaker.panels.eos.NewEOSDialog;
import ca.tonita.physics.gr.hydro.TOVBuilder;
import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.StandardTickUnitSource;

/**
 *
 * @author atonita
 */
public class TOVFamilyPanel extends javax.swing.JPanel implements ChangeListener {

    JFreeChart chart = null;
    private JawBreakerModel model;
    private TOVFamilyDataset tovFamilyDataset;
    ArrayList<JCheckBox> eosCheckBoxes = new ArrayList<JCheckBox>();

    /**
     * Creates new form TOVBuilderPanel
     */
    public TOVFamilyPanel() {
        model = new JawBreakerModel();
        setModel(model);
        // Set up GUI components.
        initComponents();
        chartPanel.setLayout(new BorderLayout());
        initChart();
    }

    public void setModel(JawBreakerModel model) {
        // Set up links to model.
        this.model = model;
        model.addEOSChangeListener(this);
        tovFamilyDataset = model.getTOVFamilies();
        if (chart != null) {
            chartPanel.removeAll();
            initChart();
        }
    }

    private void updateChart() {
        int iX = domainComboBox.getSelectedIndex();
        int iY = rangeComboBox.getSelectedIndex();
        tovFamilyDataset.setDomainVariable(iX);
        tovFamilyDataset.setRangeVariable(iY);
        String domainName = tovFamilyDataset.getDomainName();
        String rangeName = tovFamilyDataset.getRangeName();
        if (logarithmDomain.isSelected()) {
            LogarithmicAxis axis = new LogarithmicAxis("Log(" + domainName + ")");
            axis.setStandardTickUnits(new StandardTickUnitSource());
            chart.getXYPlot().setDomainAxis(axis);
        } else {
            NumberAxis axis = new NumberAxis(domainName);
            axis.setStandardTickUnits(new StandardTickUnitSource());
            chart.getXYPlot().setDomainAxis(axis);
        }
        if (logarithmRange.isSelected()) {
            LogarithmicAxis axis = new LogarithmicAxis("Log(" + rangeName + ")");
            axis.setStandardTickUnits(new StandardTickUnitSource());
            chart.getXYPlot().setRangeAxis(axis);
        } else {
            NumberAxis axis = new NumberAxis(rangeName);
            axis.setStandardTickUnits(new StandardTickUnitSource());
            chart.getXYPlot().setRangeAxis(axis);
        }
        chart.getXYPlot().datasetChanged(null);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        chartControlSplitPane = new javax.swing.JSplitPane();
        chartPanel = new javax.swing.JPanel();
        familyControlPanel = new javax.swing.JTabbedPane();
        familyCreationControlPanel = new javax.swing.JPanel();
        eosLabel = new javax.swing.JLabel();
        eosComboBox = new javax.swing.JComboBox();
        rkPanel = new javax.swing.JPanel();
        stepSizeLabel = new javax.swing.JLabel();
        stepSizeField = new javax.swing.JFormattedTextField();
        minPressureLabel = new javax.swing.JLabel();
        minPressureField = new javax.swing.JFormattedTextField();
        outputEveryLabel = new javax.swing.JLabel();
        outputEveryField = new javax.swing.JFormattedTextField();
        familyParameterPanel = new javax.swing.JPanel();
        familyMinPressureLabel = new javax.swing.JLabel();
        familyMaxPressureLabel = new javax.swing.JLabel();
        matchEOSPressureButton = new javax.swing.JButton();
        familyMinPressureField = new javax.swing.JFormattedTextField();
        familyMaxPressureField = new javax.swing.JFormattedTextField();
        familyNPointsLabel = new javax.swing.JLabel();
        familyNPointsSpinner = new javax.swing.JSpinner();
        spacingComboBox = new javax.swing.JComboBox();
        createFamilyButton = new javax.swing.JButton();
        familyExplorerControlPanel = new javax.swing.JPanel();
        rangeDomainPanel = new javax.swing.JPanel();
        rangeComboBox = new javax.swing.JComboBox();
        rangeLabel = new javax.swing.JLabel();
        domainComboBox = new javax.swing.JComboBox();
        domainLabel = new javax.swing.JLabel();
        logarithmDomain = new javax.swing.JCheckBox();
        logarithmRange = new javax.swing.JCheckBox();
        familyDisplayOuterPanel = new javax.swing.JPanel();
        familyDisplayPanel = new javax.swing.JPanel();

        setMaximumSize(new java.awt.Dimension(1224, 768));
        setMinimumSize(new java.awt.Dimension(1224, 768));
        setName(""); // NOI18N

        chartControlSplitPane.setDividerLocation(860);

        javax.swing.GroupLayout chartPanelLayout = new javax.swing.GroupLayout(chartPanel);
        chartPanel.setLayout(chartPanelLayout);
        chartPanelLayout.setHorizontalGroup(
            chartPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 859, Short.MAX_VALUE)
        );
        chartPanelLayout.setVerticalGroup(
            chartPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 842, Short.MAX_VALUE)
        );

        chartControlSplitPane.setLeftComponent(chartPanel);

        eosLabel.setText("Equation of State");

        rkPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Runge-Kutta 4 Parameters"));

        stepSizeLabel.setText("Step Size:");

        stepSizeField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        stepSizeField.setText("1.0E-3");

        minPressureLabel.setText("Minimum pressure:");

        minPressureField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        minPressureField.setText("1.5E-8");

        outputEveryLabel.setText("Save data every N steps:");

        outputEveryField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter(new java.text.DecimalFormat("###0"))));
        outputEveryField.setText("1");

        javax.swing.GroupLayout rkPanelLayout = new javax.swing.GroupLayout(rkPanel);
        rkPanel.setLayout(rkPanelLayout);
        rkPanelLayout.setHorizontalGroup(
            rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rkPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(stepSizeField)
                    .addComponent(minPressureField, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(rkPanelLayout.createSequentialGroup()
                        .addGroup(rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(stepSizeLabel)
                            .addComponent(minPressureLabel)
                            .addComponent(outputEveryLabel))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(outputEveryField))
                .addContainerGap())
        );
        rkPanelLayout.setVerticalGroup(
            rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rkPanelLayout.createSequentialGroup()
                .addComponent(stepSizeLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(stepSizeField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(minPressureLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(minPressureField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(outputEveryLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(outputEveryField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        familyParameterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Family parameters"));

        familyMinPressureLabel.setText("Minimum pressure");

        familyMaxPressureLabel.setText("Maximum pressure");

        matchEOSPressureButton.setText("Get pressure maximum from table");
        matchEOSPressureButton.setToolTipText("");
        matchEOSPressureButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                matchEOSPressureButtonActionPerformed(evt);
            }
        });

        familyMinPressureField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        familyMinPressureField.setText("1.0E-5");

        familyMaxPressureField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        familyMaxPressureField.setText("1.5E-1");

        familyNPointsLabel.setText("Number of TOVs");

        familyNPointsSpinner.setValue(50);

        spacingComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Linearly spaced", "Logarithmically spaced", "Quadraticly spaced" }));
        spacingComboBox.setSelectedIndex(1);

        javax.swing.GroupLayout familyParameterPanelLayout = new javax.swing.GroupLayout(familyParameterPanel);
        familyParameterPanel.setLayout(familyParameterPanelLayout);
        familyParameterPanelLayout.setHorizontalGroup(
            familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(familyParameterPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(familyParameterPanelLayout.createSequentialGroup()
                        .addGroup(familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, familyParameterPanelLayout.createSequentialGroup()
                                .addComponent(familyNPointsLabel)
                                .addGap(18, 18, 18)
                                .addComponent(familyNPointsSpinner, javax.swing.GroupLayout.PREFERRED_SIZE, 183, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(matchEOSPressureButton)
                            .addGroup(familyParameterPanelLayout.createSequentialGroup()
                                .addGroup(familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(familyMaxPressureLabel)
                                    .addComponent(familyMinPressureLabel))
                                .addGap(18, 18, 18)
                                .addGroup(familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addComponent(familyMinPressureField, javax.swing.GroupLayout.DEFAULT_SIZE, 169, Short.MAX_VALUE)
                                    .addComponent(familyMaxPressureField))))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(spacingComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        familyParameterPanelLayout.setVerticalGroup(
            familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(familyParameterPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(familyMinPressureLabel)
                    .addComponent(familyMinPressureField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(familyMaxPressureLabel)
                    .addComponent(familyMaxPressureField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(matchEOSPressureButton)
                .addGap(18, 18, 18)
                .addGroup(familyParameterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(familyNPointsLabel)
                    .addComponent(familyNPointsSpinner, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(spacingComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        createFamilyButton.setText("Create TOV family");
        createFamilyButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                createFamilyButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout familyCreationControlPanelLayout = new javax.swing.GroupLayout(familyCreationControlPanel);
        familyCreationControlPanel.setLayout(familyCreationControlPanelLayout);
        familyCreationControlPanelLayout.setHorizontalGroup(
            familyCreationControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(familyCreationControlPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(familyCreationControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(familyCreationControlPanelLayout.createSequentialGroup()
                        .addGroup(familyCreationControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(eosComboBox, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addGroup(familyCreationControlPanelLayout.createSequentialGroup()
                                .addComponent(eosLabel)
                                .addGap(0, 0, Short.MAX_VALUE)))
                        .addGap(7, 7, 7))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, familyCreationControlPanelLayout.createSequentialGroup()
                        .addGroup(familyCreationControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(familyParameterPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(rkPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addContainerGap())
                    .addGroup(familyCreationControlPanelLayout.createSequentialGroup()
                        .addComponent(createFamilyButton)
                        .addGap(0, 0, Short.MAX_VALUE))))
        );
        familyCreationControlPanelLayout.setVerticalGroup(
            familyCreationControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(familyCreationControlPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(eosLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(eosComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(rkPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(familyParameterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(createFamilyButton)
                .addContainerGap(228, Short.MAX_VALUE))
        );

        familyControlPanel.addTab("Family creation", familyCreationControlPanel);

        rangeDomainPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Plot properties"));

        rangeComboBox.setModel(new DefaultComboBoxModel(tovFamilyDataset.getVariableNames()));
        rangeComboBox.setSelectedIndex(2);
        rangeComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rangeComboBoxActionPerformed(evt);
            }
        });

        rangeLabel.setText("Range Variable:");

        domainComboBox.setModel(new DefaultComboBoxModel(tovFamilyDataset.getVariableNames()));
        domainComboBox.setSelectedIndex(1);
        domainComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                domainComboBoxActionPerformed(evt);
            }
        });

        domainLabel.setText("Domain Variable:");

        logarithmDomain.setSelected(true);
        logarithmDomain.setText("logarithmic");
        logarithmDomain.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                logarithmDomainActionPerformed(evt);
            }
        });

        logarithmRange.setText("logarithmic");
        logarithmRange.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                logarithmRangeActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout rangeDomainPanelLayout = new javax.swing.GroupLayout(rangeDomainPanel);
        rangeDomainPanel.setLayout(rangeDomainPanelLayout);
        rangeDomainPanelLayout.setHorizontalGroup(
            rangeDomainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rangeDomainPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(rangeDomainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(rangeComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(domainComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(rangeDomainPanelLayout.createSequentialGroup()
                        .addComponent(domainLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 106, Short.MAX_VALUE)
                        .addComponent(logarithmDomain))
                    .addGroup(rangeDomainPanelLayout.createSequentialGroup()
                        .addComponent(rangeLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(logarithmRange)))
                .addContainerGap())
        );
        rangeDomainPanelLayout.setVerticalGroup(
            rangeDomainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rangeDomainPanelLayout.createSequentialGroup()
                .addContainerGap(13, Short.MAX_VALUE)
                .addGroup(rangeDomainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(domainLabel)
                    .addComponent(logarithmDomain))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(domainComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(rangeDomainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rangeLabel)
                    .addComponent(logarithmRange))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(rangeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        familyDisplayOuterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Available families"));
        familyDisplayOuterPanel.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEADING));

        familyDisplayPanel.setLayout(new java.awt.GridBagLayout());
        familyDisplayOuterPanel.add(familyDisplayPanel);

        javax.swing.GroupLayout familyExplorerControlPanelLayout = new javax.swing.GroupLayout(familyExplorerControlPanel);
        familyExplorerControlPanel.setLayout(familyExplorerControlPanelLayout);
        familyExplorerControlPanelLayout.setHorizontalGroup(
            familyExplorerControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, familyExplorerControlPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(familyExplorerControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(familyDisplayOuterPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(rangeDomainPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        familyExplorerControlPanelLayout.setVerticalGroup(
            familyExplorerControlPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, familyExplorerControlPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(familyDisplayOuterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 529, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 85, Short.MAX_VALUE)
                .addComponent(rangeDomainPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        familyControlPanel.addTab("Family explorer", familyExplorerControlPanel);

        chartControlSplitPane.setRightComponent(familyControlPanel);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(chartControlSplitPane)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(chartControlSplitPane)
        );
    }// </editor-fold>//GEN-END:initComponents

    private void rangeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rangeComboBoxActionPerformed
        updateChart();
    }//GEN-LAST:event_rangeComboBoxActionPerformed

    private void domainComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_domainComboBoxActionPerformed
        updateChart();
    }//GEN-LAST:event_domainComboBoxActionPerformed

    private void logarithmDomainActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_logarithmDomainActionPerformed
        updateChart();
    }//GEN-LAST:event_logarithmDomainActionPerformed

    private void logarithmRangeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_logarithmRangeActionPerformed
        updateChart();
    }//GEN-LAST:event_logarithmRangeActionPerformed

    private void matchEOSPressureButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_matchEOSPressureButtonActionPerformed
        TabulatedHermite eos = model.getEos(eosComboBox.getSelectedIndex());
        double[] pressureExtrema = eos.getPressureExtrema();
        familyMaxPressureField.setText(pressureExtrema[1] + "");
    }//GEN-LAST:event_matchEOSPressureButtonActionPerformed

    private void createFamilyButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_createFamilyButtonActionPerformed
        // Get EOS and TOV Family.
        int i = eosComboBox.getSelectedIndex();
        TOVFamily family = model.getTOVFamily(i);
        family.clear();

        // Get the parameters of the family to create.
        double minFamilyPressure = Double.valueOf(familyMinPressureField.getText());
        double maxFamilyPressure = Double.valueOf(familyMaxPressureField.getText());
        int nTOVs = (Integer) familyNPointsSpinner.getValue();
        int method = spacingComboBox.getSelectedIndex();

        // Get the RK4 parameters.
        double stepSize = Double.valueOf(stepSizeField.getText());
        int outputEvery = Integer.valueOf(outputEveryField.getText());
        double minPressure = Double.valueOf(minPressureField.getText());

        if (minPressure > minFamilyPressure) {
            JOptionPane.showMessageDialog(this, "Family inconsistent", "The minimum pressure of the RK4 algorithm should be\nless than the minimum pressure of the family.", JOptionPane.WARNING_MESSAGE);
        } else {
            try {
                TOVBuilder.fillFamily(family, nTOVs, minFamilyPressure, maxFamilyPressure, method, stepSize, outputEvery, minPressure, true);
            } catch (IllegalArgumentException e) {
                JOptionPane.showMessageDialog(this, e.getMessage(), "Family creation error.", JOptionPane.ERROR_MESSAGE);
            }
        }
        updateChart();
    }//GEN-LAST:event_createFamilyButtonActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JSplitPane chartControlSplitPane;
    private javax.swing.JPanel chartPanel;
    private javax.swing.JButton createFamilyButton;
    private javax.swing.JComboBox domainComboBox;
    private javax.swing.JLabel domainLabel;
    private javax.swing.JComboBox eosComboBox;
    private javax.swing.JLabel eosLabel;
    private javax.swing.JTabbedPane familyControlPanel;
    private javax.swing.JPanel familyCreationControlPanel;
    private javax.swing.JPanel familyDisplayOuterPanel;
    private javax.swing.JPanel familyDisplayPanel;
    private javax.swing.JPanel familyExplorerControlPanel;
    private javax.swing.JFormattedTextField familyMaxPressureField;
    private javax.swing.JLabel familyMaxPressureLabel;
    private javax.swing.JFormattedTextField familyMinPressureField;
    private javax.swing.JLabel familyMinPressureLabel;
    private javax.swing.JLabel familyNPointsLabel;
    private javax.swing.JSpinner familyNPointsSpinner;
    private javax.swing.JPanel familyParameterPanel;
    private javax.swing.JCheckBox logarithmDomain;
    private javax.swing.JCheckBox logarithmRange;
    private javax.swing.JButton matchEOSPressureButton;
    private javax.swing.JFormattedTextField minPressureField;
    private javax.swing.JLabel minPressureLabel;
    private javax.swing.JFormattedTextField outputEveryField;
    private javax.swing.JLabel outputEveryLabel;
    private javax.swing.JComboBox rangeComboBox;
    private javax.swing.JPanel rangeDomainPanel;
    private javax.swing.JLabel rangeLabel;
    private javax.swing.JPanel rkPanel;
    private javax.swing.JComboBox spacingComboBox;
    private javax.swing.JFormattedTextField stepSizeField;
    private javax.swing.JLabel stepSizeLabel;
    // End of variables declaration//GEN-END:variables

    public void stateChanged(ChangeEvent e) {
        eosComboBox.setModel(new DefaultComboBoxModel(model.getEOSNames()));
        String[] eosNames = model.getEOSNames();
        // Get currently active eosModels.
        ArrayList<String> inactive = new ArrayList<String>();
        for (int i = 0; i < eosCheckBoxes.size(); i++) {
            if (!eosCheckBoxes.get(i).isSelected()) {
                inactive.add(eosCheckBoxes.get(i).getText());
            }
        }
        // Empty and refill eosCheckBoxes.
        eosCheckBoxes.clear();
        for (int i = 0; i < eosNames.length; i++) {
            JCheckBox eosBox = createEOSCheckBox(eosNames[i]);
            if (inactive.contains(eosNames[i])) {
                eosBox.setSelected(false);
            } else {
                eosBox.setSelected(true);
            }
            eosCheckBoxes.add(eosBox);
        }
        renderEOSDisplay();
    }

    /**
     * Adds all the equations of state check boxes to the eosDisplayPanel.
     */
    private void renderEOSDisplay() {
        // Empty...
        familyDisplayPanel.removeAll();

        GridBagConstraints c = new GridBagConstraints();
        c.gridx = 0;
        c.gridwidth = 1;
        c.gridheight = 1;
        c.weighty = 1.0;
        c.weightx = .1;
        c.anchor = GridBagConstraints.WEST;
        c.fill = GridBagConstraints.HORIZONTAL;

        for (int i = 0; i < eosCheckBoxes.size(); i++) {
            c.gridy = i;
            familyDisplayPanel.add(eosCheckBoxes.get(i), c);
            familyDisplayOuterPanel.validate();
        }

        chart.getXYPlot().datasetChanged(null);
    }

    /**
     * Creates an equation of state check box adding this as an action listener.
     *
     * @param id the name of the equation of state.
     * @return the created checkbox
     */
    private JCheckBox createEOSCheckBox(String id) {
        JCheckBox eosBox = new JCheckBox(id);
        eosBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                int index = eosCheckBoxes.indexOf(e.getSource());
                tovFamilyDataset.setActivated(index, eosCheckBoxes.get(index).isSelected());
                chart.getXYPlot().datasetChanged(null);
            }
        });
        eosBox.setSelected(false);
        return eosBox;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        try {
            javax.swing.UIManager.setLookAndFeel(javax.swing.UIManager.getSystemLookAndFeelClassName());
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(NewEOSDialog.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            Logger.getLogger(NewEOSDialog.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            Logger.getLogger(NewEOSDialog.class.getName()).log(Level.SEVERE, null, ex);
        } catch (UnsupportedLookAndFeelException ex) {
            Logger.getLogger(NewEOSDialog.class.getName()).log(Level.SEVERE, null, ex);
        }

        /*
         * Create and display the form
         */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                JFrame aFrame = new JFrame();
                aFrame.getContentPane().setLayout(new BorderLayout());
                aFrame.setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
                TOVFamilyPanel tovFamilyPanel = new TOVFamilyPanel();
                aFrame.add(tovFamilyPanel, BorderLayout.CENTER);
                aFrame.setVisible(true);
                aFrame.setSize(1024, 768);
                aFrame.pack();
            }
        });
    }

    private void initChart() {
        ChartPanel panel = ChartPanelCreator.createChartPanel("TOV Family", tovFamilyDataset.getDomainName(), tovFamilyDataset.getRangeName(), tovFamilyDataset);
        chart = panel.getChart();
        chartPanel.add(panel, BorderLayout.CENTER);
    }
}
