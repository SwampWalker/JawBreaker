package ca.tonita.jawbreaker.datasets;

/**
 *
 * @author atonita
 */
public class EOSBean {

    public double[][] getData() {
        return data;
    }

    public void setData(double[][] data) {
        this.data = data;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    private double[][] data;
    private String name;
}
