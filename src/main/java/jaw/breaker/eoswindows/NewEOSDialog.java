/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jaw.breaker.eoswindows;

import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.UnsupportedLookAndFeelException;
import jaw.breaker.equationsOfState.TabulatedHermite;

/**
 * Dialog for adding new equations of state.
 *
 * @author tonita
 */
public class NewEOSDialog extends javax.swing.JDialog {

    EOSStorage storage;

    /**
     * Sets the equation of state storage.
     *
     * @param storage The storage object to set.
     */
    public void setStorage(EOSStorage storage) {
        this.storage = storage;
    }

    /**
     * Creates new form NewEOSDialog with storage
     * @param parent The parent window.
     * @param modal whether this window blocks input or not
     * @param storage <code>EOSStorage</code> for holding equations of state.
     */
    public NewEOSDialog(java.awt.Frame parent, boolean modal, EOSStorage storage) {
        super(parent, modal);
        this.storage = storage;
        initComponents();
    }

    /**
     * Creates new form NewEOSDialog
     */
    public NewEOSDialog(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
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

        jTabbedPane1 = new javax.swing.JTabbedPane();
        newTableFromFilePanel1 = new jaw.breaker.eoswindows.NewTableFromFilePanel();
        newPolyPanel1 = new jaw.breaker.eoswindows.NewPolyPanel();
        jButton1 = new javax.swing.JButton();
        cancelButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setAlwaysOnTop(true);
        setMaximumSize(new java.awt.Dimension(345, 368));
        setMinimumSize(new java.awt.Dimension(345, 368));
        setPreferredSize(new java.awt.Dimension(345, 368));
        setResizable(false);

        jTabbedPane1.addTab("Tabulated EOS", newTableFromFilePanel1);
        jTabbedPane1.addTab("New Polytrope", newPolyPanel1);

        jButton1.setText("Create Equation of State table");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancelButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jTabbedPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 327, Short.MAX_VALUE)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jButton1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(cancelButton)
                .addContainerGap(30, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jTabbedPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jButton1)
                    .addComponent(cancelButton)))
        );

        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setBounds((screenSize.width-343)/2, (screenSize.height-402)/2, 343, 402);
    }// </editor-fold>//GEN-END:initComponents

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
       TabulatedHermite eos = ((TabulatedEOSGenerator) jTabbedPane1.getSelectedComponent()).getTabulatedEOS();
        if (storage != null && eos != null) {
            storage.add(eos);
            this.setVisible(false);
        }
    }//GEN-LAST:event_jButton1ActionPerformed

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed
        this.setVisible(false);
    }//GEN-LAST:event_cancelButtonActionPerformed

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
         * Create and display the dialog
         */
        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                NewEOSDialog dialog = new NewEOSDialog(new javax.swing.JFrame(), true);
                dialog.addWindowListener(new java.awt.event.WindowAdapter() {

                    @Override
                    public void windowClosing(java.awt.event.WindowEvent e) {
                        System.exit(0);
                    }
                });
                dialog.setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton cancelButton;
    private javax.swing.JButton jButton1;
    private javax.swing.JTabbedPane jTabbedPane1;
    private jaw.breaker.eoswindows.NewPolyPanel newPolyPanel1;
    private jaw.breaker.eoswindows.NewTableFromFilePanel newTableFromFilePanel1;
    // End of variables declaration//GEN-END:variables
}
