/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.panels;

import ca.tonita.jawbreaker.datasets.TOVDataset;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.models.JawBreakerModel;
import ca.tonita.physics.gr.hydro.TOVData;
import ca.tonita.jawbreaker.panels.eos.NewEOSDialog;
import ca.tonita.physics.gr.hydro.TOVBuilder;
import java.awt.BorderLayout;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JFrame;
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
public class TOVBuilderPanel extends javax.swing.JPanel implements ChangeListener {

    JFreeChart chart = null;
    private TOVData rk4TOV = new TOVData();
    private JawBreakerModel model;
    private final TOVDataset tovDataset;

    /**
     * Creates new form TOVBuilderPanel
     */
    public TOVBuilderPanel() {
        model = new JawBreakerModel();
        this.setModel(model);
        tovDataset = new TOVDataset();

        // Set up GUI components.
        initComponents();
        jPanelLeft.setLayout(new BorderLayout());
        ChartPanel panel = ChartPanelCreator.createChartPanel("TOV Model", tovDataset.getDomainName(), tovDataset.getRangeName(), tovDataset);
        chart = panel.getChart();
        jPanelLeft.add(panel, BorderLayout.CENTER);
    }
    
    public final void setModel(JawBreakerModel model) {
        // Set up links to model.
        this.model = model;
        model.addEOSChangeListener(this);
    }

    private void updateChart() {
        int iX = domainComboBox.getSelectedIndex();
        int iY = rangeComboBox.getSelectedIndex();
        tovDataset.setDomainVariable(iX);
        tovDataset.setRangeVariable(iY);
        String domainName = tovDataset.getDomainName();
        String rangeName = tovDataset.getRangeName();
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

        jSplitPane1 = new javax.swing.JSplitPane();
        jPanelLeft = new javax.swing.JPanel();
        jPanelRight = new javax.swing.JPanel();
        eosLabel = new javax.swing.JLabel();
        eosComboBox = new javax.swing.JComboBox();
        rkPanel = new javax.swing.JPanel();
        stepSizeLabel = new javax.swing.JLabel();
        stepSizeField = new javax.swing.JFormattedTextField();
        minPressureLabel = new javax.swing.JLabel();
        minPressureField = new javax.swing.JFormattedTextField();
        outputEveryLabel = new javax.swing.JLabel();
        outputEveryField = new javax.swing.JFormattedTextField();
        centralPressureLabel = new javax.swing.JLabel();
        centralPressureField = new javax.swing.JFormattedTextField();
        jPanel1 = new javax.swing.JPanel();
        massLabel = new javax.swing.JLabel();
        massField = new javax.swing.JTextField();
        radiusLabel = new javax.swing.JLabel();
        radiusField = new javax.swing.JTextField();
        rangeComboBox = new javax.swing.JComboBox();
        rangeLabel = new javax.swing.JLabel();
        domainComboBox = new javax.swing.JComboBox();
        domainLabel = new javax.swing.JLabel();
        logarithmDomain = new javax.swing.JCheckBox();
        logarithmRange = new javax.swing.JCheckBox();
        spectralPanel = new javax.swing.JPanel();
        rankSpinner = new javax.swing.JSpinner();
        jLabel1 = new javax.swing.JLabel();
        createTOVButton = new javax.swing.JButton();

        setMaximumSize(new java.awt.Dimension(1224, 768));
        setMinimumSize(new java.awt.Dimension(1224, 768));
        setName(""); // NOI18N

        jSplitPane1.setDividerLocation(860);

        javax.swing.GroupLayout jPanelLeftLayout = new javax.swing.GroupLayout(jPanelLeft);
        jPanelLeft.setLayout(jPanelLeftLayout);
        jPanelLeftLayout.setHorizontalGroup(
            jPanelLeftLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 859, Short.MAX_VALUE)
        );
        jPanelLeftLayout.setVerticalGroup(
            jPanelLeftLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 766, Short.MAX_VALUE)
        );

        jSplitPane1.setLeftComponent(jPanelLeft);

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
                        .addGap(0, 158, Short.MAX_VALUE))
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

        centralPressureLabel.setText("Desired central pressure:");

        centralPressureField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        centralPressureField.setText("7.5e-2");

        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Properties"));

        massLabel.setText("Gravitational Mass");

        massField.setEnabled(false);

        radiusLabel.setText("Radius");

        radiusField.setEnabled(false);

        rangeComboBox.setModel(new DefaultComboBoxModel(tovDataset.getVariableNames()));
        rangeComboBox.setSelectedIndex(1);
        rangeComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rangeComboBoxActionPerformed(evt);
            }
        });

        rangeLabel.setText("Range Variable:");

        domainComboBox.setModel(new DefaultComboBoxModel(tovDataset.getVariableNames()));
        domainComboBox.setSelectedIndex(0);
        domainComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                domainComboBoxActionPerformed(evt);
            }
        });

        domainLabel.setText("Domain Variable:");

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

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(rangeComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(domainComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(massLabel)
                            .addComponent(massField, javax.swing.GroupLayout.PREFERRED_SIZE, 126, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(radiusLabel)
                            .addComponent(radiusField, javax.swing.GroupLayout.DEFAULT_SIZE, 126, Short.MAX_VALUE))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(domainLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 111, Short.MAX_VALUE)
                        .addComponent(logarithmDomain))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(rangeLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(logarithmRange)))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(massLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(massField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(radiusLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(radiusField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 12, Short.MAX_VALUE)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(domainLabel)
                    .addComponent(logarithmDomain))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(domainComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rangeLabel)
                    .addComponent(logarithmRange))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(rangeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        spectralPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Spectral Properties"));

        rankSpinner.setModel(new javax.swing.SpinnerNumberModel(20, 1, 255, 1));

        jLabel1.setText("Number of abscissas");

        javax.swing.GroupLayout spectralPanelLayout = new javax.swing.GroupLayout(spectralPanel);
        spectralPanel.setLayout(spectralPanelLayout);
        spectralPanelLayout.setHorizontalGroup(
            spectralPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, spectralPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(rankSpinner, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        spectralPanelLayout.setVerticalGroup(
            spectralPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spectralPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(spectralPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rankSpinner, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel1))
                .addContainerGap(45, Short.MAX_VALUE))
        );

        createTOVButton.setText("Make Model");
        createTOVButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                createTOVButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanelRightLayout = new javax.swing.GroupLayout(jPanelRight);
        jPanelRight.setLayout(jPanelRightLayout);
        jPanelRightLayout.setHorizontalGroup(
            jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanelRightLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spectralPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPanel1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(eosComboBox, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(rkPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(centralPressureField)
                    .addGroup(jPanelRightLayout.createSequentialGroup()
                        .addGroup(jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(eosLabel)
                            .addComponent(centralPressureLabel)
                            .addComponent(createTOVButton))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanelRightLayout.setVerticalGroup(
            jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanelRightLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(eosLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(eosComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(12, 12, 12)
                .addComponent(centralPressureLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(centralPressureField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(24, 24, 24)
                .addComponent(rkPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(spectralPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(createTOVButton)
                .addGap(24, 24, 24)
                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        jSplitPane1.setRightComponent(jPanelRight);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1)
        );
    }// </editor-fold>//GEN-END:initComponents

    private void createTOVButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_createTOVButtonActionPerformed
        double centralPressure = Double.valueOf(centralPressureField.getText());
        TabulatedHermite eos = model.getEos(eosComboBox.getSelectedIndex());
        double stepSize = Double.valueOf(stepSizeField.getText());
        int outputEvery = Integer.valueOf(outputEveryField.getText());
        double minPressure = Double.valueOf(minPressureField.getText());
        TOVBuilder.evolve(rk4TOV, eos, centralPressure, stepSize, outputEvery, minPressure);
        rk4TOV.computeSecondaries(eos);
        if (tovDataset.getSeriesCount() == 0) {
            tovDataset.add(0, rk4TOV);
        }
        massField.setText(String.format("%.6f", rk4TOV.getConservedMass()));
        radiusField.setText(String.format("%.6f",rk4TOV.getRadius()));
        chart.getXYPlot().datasetChanged(null);
    }//GEN-LAST:event_createTOVButtonActionPerformed

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
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JFormattedTextField centralPressureField;
    private javax.swing.JLabel centralPressureLabel;
    private javax.swing.JButton createTOVButton;
    private javax.swing.JComboBox domainComboBox;
    private javax.swing.JLabel domainLabel;
    private javax.swing.JComboBox eosComboBox;
    private javax.swing.JLabel eosLabel;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanelLeft;
    private javax.swing.JPanel jPanelRight;
    private javax.swing.JSplitPane jSplitPane1;
    private javax.swing.JCheckBox logarithmDomain;
    private javax.swing.JCheckBox logarithmRange;
    private javax.swing.JTextField massField;
    private javax.swing.JLabel massLabel;
    private javax.swing.JFormattedTextField minPressureField;
    private javax.swing.JLabel minPressureLabel;
    private javax.swing.JFormattedTextField outputEveryField;
    private javax.swing.JLabel outputEveryLabel;
    private javax.swing.JTextField radiusField;
    private javax.swing.JLabel radiusLabel;
    private javax.swing.JComboBox rangeComboBox;
    private javax.swing.JLabel rangeLabel;
    private javax.swing.JSpinner rankSpinner;
    private javax.swing.JPanel rkPanel;
    private javax.swing.JPanel spectralPanel;
    private javax.swing.JFormattedTextField stepSizeField;
    private javax.swing.JLabel stepSizeLabel;
    // End of variables declaration//GEN-END:variables

    public void stateChanged(ChangeEvent e) {
        eosComboBox.setModel(new DefaultComboBoxModel(model.getEOSNames()));
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
                TOVBuilderPanel tovBuilderPanel = new TOVBuilderPanel();
                aFrame.add(tovBuilderPanel, BorderLayout.CENTER);
                aFrame.setVisible(true);
                aFrame.setSize(1024, 768);
                aFrame.pack();
            }
        });
    }
}
