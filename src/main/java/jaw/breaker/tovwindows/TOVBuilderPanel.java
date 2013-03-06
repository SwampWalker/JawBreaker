/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jaw.breaker.tovwindows;

import java.awt.BorderLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import jaw.breaker.datasets.TOVDataset;
import jaw.breaker.eoswindows.ChartPanelCreator;
import jaw.breaker.eoswindows.EOSStorage;
import jaw.breaker.models.TOVData;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

/**
 *
 * @author atonita
 */
public class TOVBuilderPanel extends javax.swing.JPanel implements ChangeListener {
    JFreeChart chart = null;
    private EOSStorage eosStorage;
    private TOVData tov = new TOVData();

    /**
     * Creates new form TOVBuilderPanel
     */
    public TOVBuilderPanel() {
        TOVDataset tovData = new TOVDataset();
        tovData.setTOV(tov);
        initComponents();
        jPanelLeft.setLayout(new BorderLayout());
        ChartPanel panel = ChartPanelCreator.createChartPanel("TOV Model", tovData.getDomainName(), tovData.getRangeName(), tovData);
        chart = panel.getChart();
        jPanelLeft.add(panel, BorderLayout.CENTER);
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
        maxStepsLabel = new javax.swing.JLabel();
        maxStepsField = new javax.swing.JFormattedTextField();
        maxRadiusLabel = new javax.swing.JLabel();
        maxRadiusField = new javax.swing.JFormattedTextField();
        minPressureLabel = new javax.swing.JLabel();
        minPressureField = new javax.swing.JFormattedTextField();
        centralPressureLabel = new javax.swing.JLabel();
        centralPressureField = new javax.swing.JFormattedTextField();
        jButton1 = new javax.swing.JButton();

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

        maxStepsLabel.setText("Maximum number of steps:");

        maxStepsField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter(java.text.NumberFormat.getIntegerInstance())));
        maxStepsField.setText("10000");

        maxRadiusLabel.setText("Maximum radius:");

        maxRadiusField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        maxRadiusField.setText("1.5E1");

        minPressureLabel.setText("Minimum pressure:");

        minPressureField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        minPressureField.setText("1.5E-8");

        javax.swing.GroupLayout rkPanelLayout = new javax.swing.GroupLayout(rkPanel);
        rkPanel.setLayout(rkPanelLayout);
        rkPanelLayout.setHorizontalGroup(
            rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rkPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(stepSizeField)
                    .addComponent(maxStepsField)
                    .addComponent(maxRadiusField, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(rkPanelLayout.createSequentialGroup()
                        .addGroup(rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(stepSizeLabel)
                            .addComponent(maxStepsLabel)
                            .addComponent(maxRadiusLabel)
                            .addComponent(minPressureLabel))
                        .addGap(0, 141, Short.MAX_VALUE))
                    .addComponent(minPressureField, javax.swing.GroupLayout.Alignment.TRAILING))
                .addContainerGap())
        );
        rkPanelLayout.setVerticalGroup(
            rkPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rkPanelLayout.createSequentialGroup()
                .addComponent(stepSizeLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(stepSizeField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(13, 13, 13)
                .addComponent(maxStepsLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(maxStepsField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(12, 12, 12)
                .addComponent(maxRadiusLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(maxRadiusField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(minPressureLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(minPressureField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        centralPressureLabel.setText("Desired central pressure:");

        centralPressureField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        centralPressureField.setText("7.5e-2");

        jButton1.setText("Make Model");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanelRightLayout = new javax.swing.GroupLayout(jPanelRight);
        jPanelRight.setLayout(jPanelRightLayout);
        jPanelRightLayout.setHorizontalGroup(
            jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanelRightLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(eosComboBox, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(rkPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(centralPressureField)
                    .addGroup(jPanelRightLayout.createSequentialGroup()
                        .addGroup(jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(eosLabel)
                            .addComponent(centralPressureLabel)
                            .addComponent(jButton1))
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
                .addComponent(rkPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 253, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton1)
                .addContainerGap(342, Short.MAX_VALUE))
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

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        double[] y0 = new double[] {Double.valueOf(centralPressureField.getText()), 0.0, 0.0};
        
    }//GEN-LAST:event_jButton1ActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JFormattedTextField centralPressureField;
    private javax.swing.JLabel centralPressureLabel;
    private javax.swing.JComboBox eosComboBox;
    private javax.swing.JLabel eosLabel;
    private javax.swing.JButton jButton1;
    private javax.swing.JPanel jPanelLeft;
    private javax.swing.JPanel jPanelRight;
    private javax.swing.JSplitPane jSplitPane1;
    private javax.swing.JFormattedTextField maxRadiusField;
    private javax.swing.JLabel maxRadiusLabel;
    private javax.swing.JFormattedTextField maxStepsField;
    private javax.swing.JLabel maxStepsLabel;
    private javax.swing.JFormattedTextField minPressureField;
    private javax.swing.JLabel minPressureLabel;
    private javax.swing.JPanel rkPanel;
    private javax.swing.JFormattedTextField stepSizeField;
    private javax.swing.JLabel stepSizeLabel;
    // End of variables declaration//GEN-END:variables

    public void setEOSStorage(EOSStorage eosStorage) {
        this.eosStorage = eosStorage;
        eosStorage.addChangeListener(this);
    }

    public void stateChanged(ChangeEvent e) {
        eosComboBox.setModel(new DefaultComboBoxModel(eosStorage.getNames()));
    }
}
