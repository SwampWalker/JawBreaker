/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.tonita.jawbreaker.panels;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import ca.tonita.jawbreaker.datasets.EOSDataset;
import ca.tonita.jawbreaker.equationsOfState.TabulatedHermite;
import ca.tonita.jawbreaker.models.JawBreakerModel;
import ca.tonita.jawbreaker.panels.eos.NewEOSDialog;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.StandardTickUnitSource;

/**
 *
 * @author tonita
 */
public class EOSPanel extends javax.swing.JPanel implements ChangeListener {

    private EOSDataset eosDataset = null;
    private JawBreakerModel model = null;
    private ArrayList<ChangeListener> changeListeners = new ArrayList<ChangeListener>();
    private JFreeChart chart;

    /**
     * Creates new form EOSPanel
     */
    public EOSPanel() {
        eosDataset = new EOSDataset();
        model = new JawBreakerModel();
        this.setModel(model);
        
        initComponents();
        jPanelLeft.setLayout(new BorderLayout());
        ChartPanel panel = ChartPanelCreator.createChartPanel("Equations of State", eosDataset.getDomainName(), eosDataset.getRangeName(), eosDataset);
        chart = panel.getChart();
        jPanelLeft.add(panel, BorderLayout.CENTER);
        plotPropertyChangeHandler(null);
    }
    
    public void setModel(JawBreakerModel model) {
        this.model = model;
        model.addEOSChangeListener(this);
        model.setEOSDataset(eosDataset);
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
                EOSPanel eosPanel = new EOSPanel();
                aFrame.add(eosPanel, BorderLayout.CENTER);
                aFrame.setVisible(true);
                aFrame.setSize(1024, 768);
                aFrame.pack();
            }
        });
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
        jPanelRight = new javax.swing.JPanel();
        newEOSButton = new javax.swing.JButton();
        eosDisplayOuterPanel = new javax.swing.JPanel();
        eosDisplayPanel = new javax.swing.JPanel();
        jPanel1 = new javax.swing.JPanel();
        domainComboBox = new javax.swing.JComboBox();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        rangeComboBox = new javax.swing.JComboBox();
        logarithmDomain = new javax.swing.JCheckBox();
        logarithmRange = new javax.swing.JCheckBox();
        deleteUnselectedButton = new javax.swing.JButton();
        jPanelLeft = new javax.swing.JPanel();

        setMaximumSize(new java.awt.Dimension(1224, 768));
        setMinimumSize(new java.awt.Dimension(1224, 768));
        setPreferredSize(new java.awt.Dimension(1224, 768));

        jSplitPane1.setDividerLocation(860);

        newEOSButton.setText("New Equation of State");
        newEOSButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                newEOSButtonActionPerformed(evt);
            }
        });

        eosDisplayOuterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Equations of state"));
        eosDisplayOuterPanel.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEADING));

        eosDisplayPanel.setLayout(new java.awt.GridBagLayout());
        eosDisplayOuterPanel.add(eosDisplayPanel);

        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Plot Properties"));

        domainComboBox.setModel(new DefaultComboBoxModel(eosDataset.getVariableNames()));
        domainComboBox.setSelectedIndex(1);
        domainComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                plotPropertyChangeHandler(evt);
            }
        });

        jLabel1.setText("Domain (x) variable:");

        jLabel2.setText("Range (y) variable:");

        rangeComboBox.setModel(new DefaultComboBoxModel(eosDataset.getVariableNames()));
        rangeComboBox.setSelectedIndex(2);
        rangeComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                plotPropertyChangeHandler(evt);
            }
        });

        logarithmDomain.setSelected(true);
        logarithmDomain.setText("Log scale in domain");
        logarithmDomain.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                plotPropertyChangeHandler(evt);
            }
        });

        logarithmRange.setSelected(true);
        logarithmRange.setText("Log scale in range");
        logarithmRange.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                plotPropertyChangeHandler(evt);
            }
        });

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(domainComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(rangeComboBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(logarithmRange)
                            .addComponent(logarithmDomain)
                            .addComponent(jLabel1)
                            .addComponent(jLabel2))
                        .addGap(0, 157, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(domainComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jLabel2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(rangeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(logarithmDomain)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(logarithmRange)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        deleteUnselectedButton.setText("Delete unselected");
        deleteUnselectedButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                deleteUnselectedButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanelRightLayout = new javax.swing.GroupLayout(jPanelRight);
        jPanelRight.setLayout(jPanelRightLayout);
        jPanelRightLayout.setHorizontalGroup(
            jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanelRightLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jPanel1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(newEOSButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(eosDisplayOuterPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 338, Short.MAX_VALUE)
                    .addGroup(jPanelRightLayout.createSequentialGroup()
                        .addComponent(deleteUnselectedButton)
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanelRightLayout.setVerticalGroup(
            jPanelRightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanelRightLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(newEOSButton)
                .addGap(18, 18, 18)
                .addComponent(eosDisplayOuterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 436, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(deleteUnselectedButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 19, Short.MAX_VALUE)
                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        jSplitPane1.setRightComponent(jPanelRight);

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

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 1224, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1)
        );
    }// </editor-fold>//GEN-END:initComponents
    NewEOSDialog newEOSDialog;
    ArrayList<TabulatedHermite> eosStorage = new ArrayList<TabulatedHermite>();
    ArrayList<JCheckBox> eosCheckBoxes = new ArrayList<JCheckBox>();

    private void newEOSButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_newEOSButtonActionPerformed
        if (newEOSDialog == null) {
            newEOSDialog = new NewEOSDialog(null, false, model);
        }
        newEOSDialog.setVisible(true);
    }//GEN-LAST:event_newEOSButtonActionPerformed

    private void plotPropertyChangeHandler(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_plotPropertyChangeHandler
        int iX = domainComboBox.getSelectedIndex();
        int iY = rangeComboBox.getSelectedIndex();
        eosDataset.setDomainVariable(iX);
        eosDataset.setRangeVariable(iY);
        String domainName = eosDataset.getDomainName();
        String rangeName = eosDataset.getRangeName();
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
    }//GEN-LAST:event_plotPropertyChangeHandler

    private void deleteUnselectedButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_deleteUnselectedButtonActionPerformed
        for (int i = eosCheckBoxes.size() - 1; i >= 0; i--) {
            if (!eosCheckBoxes.get(i).isSelected()) {
                eosCheckBoxes.remove(i);
                model.removeEOS(i);
            }
        }
    }//GEN-LAST:event_deleteUnselectedButtonActionPerformed
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton deleteUnselectedButton;
    private javax.swing.JComboBox domainComboBox;
    private javax.swing.JPanel eosDisplayOuterPanel;
    private javax.swing.JPanel eosDisplayPanel;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanelLeft;
    private javax.swing.JPanel jPanelRight;
    private javax.swing.JSplitPane jSplitPane1;
    private javax.swing.JCheckBox logarithmDomain;
    private javax.swing.JCheckBox logarithmRange;
    private javax.swing.JButton newEOSButton;
    private javax.swing.JComboBox rangeComboBox;
    // End of variables declaration//GEN-END:variables

    /**
     * Adds all the equations of state check boxes to the eosDisplayPanel.
     */
    private void renderEOSDisplay() {
        // Empty...
        eosDisplayPanel.removeAll();

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
            eosDisplayPanel.add(eosCheckBoxes.get(i), c);
            eosDisplayOuterPanel.validate();
        }

        chart.getXYPlot().datasetChanged(null);
    }

    /**
     * Receives notifications of changes to the equations of state stored by the
     * data model and updates the GUI accordingly.
     * @param e unused
     */
    public void stateChanged(ChangeEvent e) {
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
                eosDataset.setActivated(index, eosCheckBoxes.get(index).isSelected());
                chart.getXYPlot().datasetChanged(null);
            }
        });
        eosBox.setSelected(true);
        return eosBox;
    }
}
