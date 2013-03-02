/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jaw.breaker.eoswindows;

import java.text.DecimalFormat;
import java.text.ParseException;
import java.util.logging.Level;
import java.util.logging.Logger;
import jaw.breaker.equationsOfState.Polytrope;
import jaw.breaker.equationsOfState.TabulatedHermite;

/**
 *
 * @author atonita
 */
public class NewPolyPanel extends javax.swing.JPanel implements TabulatedEOSGenerator {

    /**
     * Creates new form NewPolyPanel.
     */
    public NewPolyPanel() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        gammaField = new javax.swing.JFormattedTextField();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        particleMassField = new javax.swing.JFormattedTextField();
        nPointsField = new javax.swing.JFormattedTextField();
        jLabel7 = new javax.swing.JLabel();
        logNminField = new javax.swing.JFormattedTextField();
        logNmaxField = new javax.swing.JFormattedTextField();
        jLabel8 = new javax.swing.JLabel();
        kappaField = new javax.swing.JFormattedTextField();

        jLabel1.setText("Polytropic Kappa:");

        jLabel2.setText("Adiabatic index (Γ):");

        gammaField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        gammaField.setText("2.000E0");

        jLabel3.setFont(new java.awt.Font("Tahoma", 1, 18)); // NOI18N
        jLabel3.setText("Physical properties:");

        jLabel4.setFont(new java.awt.Font("Tahoma", 1, 18)); // NOI18N
        jLabel4.setText("Table properties:");

        jLabel5.setText("Number of samples:");

        jLabel6.setText("Particle mass:");

        particleMassField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        particleMassField.setText("1.000E0");

        nPointsField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter(java.text.NumberFormat.getIntegerInstance())));
        nPointsField.setText("100");

        jLabel7.setText("Minimum log density:");
        jLabel7.setToolTipText("Base 10 logarithm of particle number density, at the initial point of the table.");

        logNminField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        logNminField.setText("-3.00E0");

        logNmaxField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        logNmaxField.setText("-1.00E0");

        jLabel8.setText("Maximum log density:");
        jLabel8.setToolTipText("Base 10 logarithm of particle number density, at the final point of the table.");

        kappaField.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(
            new javax.swing.text.NumberFormatter(
                new java.text.DecimalFormat("0.00#############E0#")
            )
        ));
        kappaField.setText("1.000E2");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(jLabel8)
                        .addGap(52, 52, 52)
                        .addComponent(logNmaxField))
                    .addComponent(jLabel3)
                    .addComponent(jLabel4)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(jLabel6, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addGap(67, 67, 67)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(gammaField)
                            .addComponent(particleMassField)))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel5)
                                    .addComponent(jLabel7))
                                .addGap(56, 56, 56))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jLabel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)))
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(nPointsField, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(logNminField, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 141, Short.MAX_VALUE)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                                .addGap(6, 6, 6)
                                .addComponent(kappaField)))))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel3)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(kappaField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(gammaField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel6)
                    .addComponent(particleMassField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(21, 21, 21)
                .addComponent(jLabel4)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel5)
                    .addComponent(nPointsField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel7)
                    .addComponent(logNminField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(logNmaxField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel8))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JFormattedTextField gammaField;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JFormattedTextField kappaField;
    private javax.swing.JFormattedTextField logNmaxField;
    private javax.swing.JFormattedTextField logNminField;
    private javax.swing.JFormattedTextField nPointsField;
    private javax.swing.JFormattedTextField particleMassField;
    // End of variables declaration//GEN-END:variables

    public TabulatedHermite getTabulatedEOS() {
        int nPoints = Integer.valueOf(nPointsField.getText());
        double logNMin = Double.valueOf(logNminField.getText());
        double logNMax = Double.valueOf(logNmaxField.getText());
        double particleMass = Double.valueOf(particleMassField.getText());
        double kappa = Double.valueOf(kappaField.getText());
        double gamma = Double.valueOf(gammaField.getText());
        
        Polytrope poly = new Polytrope(kappa, gamma, particleMass);
        double[] logn = new double[nPoints];
        double[] logp = new double[nPoints];
        double[] energyPerParticle = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            logn[i] = logNMin + i*(logNMax-logNMin)/(nPoints - 1);
            double n = Math.pow(10, logn[i]);
            logp[i] = Math.log10(poly.pressure(n));
            energyPerParticle[i] = poly.energyPerParticle(n);
        }
        TabulatedHermite eos = new TabulatedHermite(logn, logp, energyPerParticle, particleMass);
        eos.setIdentifier("Polytrope: Γ = " + gamma + ", k = " + kappa);
        return eos;
    }
}
