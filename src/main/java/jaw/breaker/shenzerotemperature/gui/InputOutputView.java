package jaw.breaker.shenzerotemperature.gui;

import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextPane;
import javax.swing.border.EmptyBorder;
import jaw.breaker.shenzerotemperature.drivers.InputOutput;

/**
 * The gui part for an InputOutput driver.
 * @author atonita
 *
 */
@SuppressWarnings("serial")
public class InputOutputView extends JFrame {

	private JPanel contentPane;
	private Action action = null;

	/**
	 * Create the frame.
	 * @param payload The payload to run when files have been selected.
	 */
	public InputOutputView(InputOutput payload, String btnLabel, String helpText) {
		action = new InputOutputController(this, payload);
		
		setTitle("Shen Zero Temperature Lepton Contributor");
		setResizable(false);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 503, 313);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(null);
		
		JButton btnAddsLeptonsTo = new JButton(btnLabel);
		btnAddsLeptonsTo.setAction(action);
		btnAddsLeptonsTo.setText(btnLabel);
		btnAddsLeptonsTo.setBounds(43, 29, 414, 43);
		contentPane.add(btnAddsLeptonsTo);
		
		JTextPane txtpnPushingThatButton = new JTextPane();
		txtpnPushingThatButton.setEditable(false);
		txtpnPushingThatButton.setText(helpText);
		txtpnPushingThatButton.setBounds(43, 113, 414, 120);
		contentPane.add(txtpnPushingThatButton);
	}
}
