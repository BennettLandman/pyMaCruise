// Derived from Sun source code:  http://java.sun.com/docs/books/tutorial/uiswing/examples/components/TextSamplerDemoProject/src/components/TextSamplerDemo.java
/*
 * Copyright (c) 1995 - 2008 Sun Microsystems, Inc.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *   - Neither the name of Sun Microsystems nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */ 

package license;
import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Stack;

import javax.swing.JButton;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

import edu.jhu.ece.iacl.jist.pipeline.gui.resources.PlaceHolder;
import edu.jhu.ece.iacl.jist.utility.JistLogger;

public class JISTLicenseAbout extends JFrame implements ActionListener {

	private static String releaseID = "@JIST-RELEASE-ID@";
	private static String releaseBuildDate = "@JIST-BUILD-DATE@";
	private JTextField url;
	private JTextField buildIDField;

	private JButton loadButton;

	private JButton backButton;

	private Stack urlStack = new Stack();

	private JEditorPane editorPane;

	public JISTLicenseAbout() {		
		setTitle("Java Image Sceience Toolkit - Release: " + releaseID);
		setSize(600, 400);
		setLayout(new BorderLayout());
		
		
		
		url = new JTextField(30);
		
		url.setEditable(false);
		//Create an editor pane.
		editorPane = createEditorPane();		
		JScrollPane editorScrollPane = new JScrollPane(editorPane);
		editorScrollPane.setVerticalScrollBarPolicy(
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

		editorScrollPane.setPreferredSize(new Dimension(250, 145));
		editorScrollPane.setMinimumSize(new Dimension(10, 10));
		add(editorScrollPane);



		loadButton = new JButton("Reload");
		loadButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent event) {
				try { // remember URL for back button
					urlStack.push(url.getText());

					editorPane.setPage(url.getText());
				} catch (IOException e) {
					editorPane.setText("Error: " + e);
				}
			}
		});

		// set up back button and button action

		backButton = new JButton("Back");
		backButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent event) {
				if (urlStack.size() <= 1)
					return;
				try { // get URL from back button
					urlStack.pop();
					// show URL in text field
					String urlString = (String) urlStack.peek();
					url.setText(urlString);

					editorPane.setPage(urlString);
				} catch (IOException e) {
					editorPane.setText("Error: " + e);
				}
			}
		});

		Container contentPane = getContentPane();
		contentPane.add(new JScrollPane(editorPane), "Center");

		// put all control components in a panel

		JPanel panel = new JPanel();
		panel.add(new JLabel("URL"));
		panel.add(url);
		panel.add(loadButton);
		panel.add(backButton);
		contentPane.add(panel, "South");
		
		buildIDField = new JTextField(10);
		buildIDField.setEditable(false);
		buildIDField.setText(releaseID);
		JPanel uppanel = new JPanel();
		uppanel .add(new JLabel("JIST Release:"));
		uppanel .add(buildIDField);
		contentPane.add(uppanel, "North");
		
		setIconImage(Toolkit.getDefaultToolkit().getImage(
				PlaceHolder.class.getResource("zoom.gif"))); /**SET ICON HERE**/
	}
	private JEditorPane createEditorPane() {	
		editorPane = new JEditorPane();
		editorPane.setEditable(false);
		editorPane.addHyperlinkListener(new HyperlinkListener() {
			public void hyperlinkUpdate(HyperlinkEvent event) {
				if (event.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
					try { // remember URL for back button
						urlStack.push(event.getURL().toString());
						// show URL in text field
						url.setText(event.getURL().toString());

						editorPane.setPage(event.getURL());
					} catch (IOException e) {
						editorPane.setText("Error: " + e);
					}
				}
			}
		});
		URL firstURL = getClass().getResource("JIST-License.html");
		
			try {
				if(firstURL==null)
				firstURL=new URL("http://yahoo.com");
			} catch (MalformedURLException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		this.url.setText(firstURL.toString());
		urlStack.push(this.url.getText());


		if (firstURL != null) {
			try {
				editorPane.setPage(firstURL);
			} catch (Exception e) {
				System.err.println(getClass().getCanonicalName()+"Attempted to read a bad URL: " + firstURL);
			}
		} else {
			System.err.println(getClass().getCanonicalName()+"Couldn't find file: TextSamplerDemoHelp.html");
		}

		return editorPane;
	}
	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		JistLogger.logOutput(JistLogger.INFO, getClass().getCanonicalName()+"\t"+arg0);

	}

	public static void main(String[] args) {
		JFrame frame = new JISTLicenseAbout();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.show();
	}
	/**
	 * Create the GUI and show it.  For thread safety,
	 * this method should be invoked from the
	 * event dispatch thread.
	 */
	/* private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("JIST License");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Add content to the window.
        frame.add(new JISTLicenseAbout());

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }*/
	public static void showDialog() {
		// TODO Auto-generated method stub
		JFrame frame = new JISTLicenseAbout();
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.show();
	}
	public static String getReleaseID() {
		if(releaseID.contains("@"))
			return "Dev";
		else 
			return releaseID;
	}
	public static String getReleaseBuildDate() {
		if(releaseBuildDate.contains("@"))
			return "private build";
		else 
			return releaseBuildDate;
	}

	/*  public static void main(String[] args) {
        //Schedule a job for the event dispatching thread:
        //creating and showing this application's GUI.
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                 //Turn off metal's use of bold fonts
		UIManager.put("swing.boldMetal", Boolean.FALSE);
		createAndShowGUI();
            }
        });
    }*/

}
