package ca.tonita.jawbreaker.shenzerotemperature.gui;

import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import javax.swing.AbstractAction;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import ca.tonita.jawbreaker.shenzerotemperature.drivers.InputOutput;

/**
 * Action event whose method is called when the button is pushed.
 * Selects files and passes them on to the InputOutput payload.
 * @author atonita
 *
 */
@SuppressWarnings("serial")
public class InputOutputController extends AbstractAction {
	/**
	 * Where the actual work will be performed.
	 */
	private InputOutput payload;
	
	/**
	 * The parent frame.
	 */
	private final JFrame display;
	
	/**
	 * Constructor.
	 * @param display The parent JFrame.
	 * @param payload
	 */
	public InputOutputController(JFrame display, InputOutput payload) {
		this.display = display;
		this.payload = payload;
		putValue(NAME, "SwingAction");
		putValue(SHORT_DESCRIPTION, "Some short description");
	}
	public void actionPerformed(ActionEvent e) {
		boolean validFiles = false;
		
		// Get equation of state file to open.
		JFileChooser fc = new JFileChooser();
		int fileState = fc.showOpenDialog(display);
		File eosFile = null;
		File outFile = null;
		if (fileState == JFileChooser.APPROVE_OPTION) {
			eosFile = fc.getSelectedFile();
			validFiles = eosFile.canRead();
			if (!validFiles) {
				JOptionPane.showMessageDialog(
						display, "Could not read file.", 
						"Read error", JOptionPane.ERROR_MESSAGE);
			} else {
				// Opening up two file selection dialogs in quick succession is confusing, so...
				JOptionPane.showMessageDialog(
						display, "File opened successfully.", 
						"Notification", JOptionPane.INFORMATION_MESSAGE);
			}
		}
		
		// If the file was good.
		if (validFiles) {
			validFiles = false;
			// Get save file
			fileState = fc.showSaveDialog(display);
			if (fileState == JFileChooser.APPROVE_OPTION) {
				outFile = fc.getSelectedFile();
				
				try {
					// See if we can create the file.
					validFiles = outFile.createNewFile();
					if (!validFiles) {
						// File exists.
						int choice = JOptionPane.showConfirmDialog(
								display, "File exists. Do you want to overwrite the file?",
								"File exists", JOptionPane.YES_NO_OPTION);
						if (choice == JOptionPane.YES_OPTION) {
							validFiles = outFile.canWrite();
							if (!validFiles) {
								JOptionPane.showMessageDialog(
										display, "Cannot write to file.", 
										"Write error", JOptionPane.ERROR_MESSAGE);
							}
						}
					}
				} catch (IOException e1) {
					JOptionPane.showMessageDialog(
							display, "I/O error.", 
							"Write error", JOptionPane.ERROR_MESSAGE);
				} catch (SecurityException e2) {
					JOptionPane.showMessageDialog(
							display, "Incorrect permissions.", 
							"Write error", JOptionPane.ERROR_MESSAGE);
				}
			}
		} // end if file opened
		
		if (validFiles) {
			payload.readAndWrite(eosFile, outFile);
		}
		
	}// end actionPerformed.
}