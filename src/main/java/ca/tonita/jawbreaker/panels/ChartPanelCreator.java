package ca.tonita.jawbreaker.panels;

import java.awt.Color;
import javax.swing.JPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.AbstractXYDataset;

/**
 *
 * @author atonita
 */
public class ChartPanelCreator {
    
    /**
     * Creates a chart panel for plotting datasets.
     *
     * @return A panel.
     */
    public static ChartPanel createChartPanel(String title, String domainLabel, String rangeLabel, AbstractXYDataset ds) {
        JFreeChart chart = ChartFactory.createScatterPlot(title, domainLabel, rangeLabel, ds, PlotOrientation.VERTICAL, true, false, false);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setNoDataMessage("NO DATA");
        plot.setDomainZeroBaselineVisible(true);
        plot.setRangeZeroBaselineVisible(true);

        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        renderer.setSeriesOutlinePaint(0, Color.black);
        renderer.setUseOutlinePaint(true);
        ChartPanel chartPanel = new ChartPanel(chart);

        chartPanel.setDomainZoomable(true);
        chartPanel.setRangeZoomable(true);
        chartPanel.setMouseWheelEnabled(true);
        return chartPanel;
    }
    
}
