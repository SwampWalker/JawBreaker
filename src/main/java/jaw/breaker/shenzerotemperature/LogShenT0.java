package jaw.breaker.shenzerotemperature;

import java.io.File;
import javax.swing.UIManager;
import jaw.breaker.shenzerotemperature.drivers.InputOutput;
import jaw.breaker.shenzerotemperature.drivers.Logarithms;
import jaw.breaker.shenzerotemperature.gui.InputOutputView;

public class LogShenT0 {

    public static void usage() {
        System.out.println("Usage: LogShenT0 [table.t00 output]");
        System.out.println("\tOpens Shen equation of state table table.t00 with leptons, "
                + "takes logarithms of baryon density and pressure and saves this to output.\n");
        System.out.println("\tCalling with no command line arguments will open the gui.");
    }
    static String btnLabel = "Take logarithm of zero-temperature equation of state";
    static String helpText = "Pushing that button will open two dialog boxes successively. "
            +"The first dialog box will ask for the zero temperature equation of state file, "
            +"as output from ShenZeroTemperature (with leptons added). The logarithm of some"
            +" quantities will be taken and saved to the second file.";

    /**
     * @param args
     */
    public static void main(String[] args) {
        InputOutput payload = new Logarithms();
        String file1 = "D:\\elasticity\\data\\eos\\eos3.leptons.yp0.t00";
        String file2 = "D:\\elasticity\\data\\eos\\eos3.log.leptons.yp0.t00";

        String[] files = {file1, file2};

        if (files.length == 2) {
            // Then don't make the gui
            File input = new File(files[0]);
            File output = new File(files[1]);
            payload.readAndWrite(input, output);
        } else if (args.length == 0) {
            // Make the gui.
            try {
                UIManager.setLookAndFeel(
                        UIManager.getSystemLookAndFeelClassName());
                InputOutputView gui = new InputOutputView(payload, btnLabel, helpText);
                gui.setVisible(true);
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            usage();
        }
    }
}
