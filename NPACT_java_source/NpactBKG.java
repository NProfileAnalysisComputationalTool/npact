import javax.swing.*;
import java.awt.*;


public class NpactBKG extends JPanel {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public void paintComponent(Graphics g)
	{
		ImageIcon icon = new ImageIcon();
		java.net.URL imageURL = NpactBKG.class.getResource("img/splash.png");
		if (imageURL != null) {
		    icon = new ImageIcon(imageURL);
		}
		
		//Image img = new ImageIcon(System.getProperty("user.dir") + "/externals/img/splash.png").getImage();
        Image img = icon.getImage();
		Dimension size = new Dimension(img.getWidth(null), img.getHeight(null));
        setPreferredSize(size);
        setMinimumSize(size);
        setMaximumSize(size);
        setSize(size);
        setLayout(null);
        g.drawImage(img, 0, 0, null);
	}
}
