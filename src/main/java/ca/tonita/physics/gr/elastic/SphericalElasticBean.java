package ca.tonita.physics.gr.elastic;

/**
 *
 * @author atonita
 */
public class SphericalElasticBean {
    private double pressure;
    private double numberDensity;
    private double massPotential;
    private double dpressure;
    private double dnumberDensity;
    private double dMassPotential;
    private double energyPerParticle;
    private double denergyPerParticle;
    private double lameLambda;
    private double shearModulus;

    public double getEnergyPerParticle() {
        return energyPerParticle;
    }

    public void setEnergyPerParticle(double energyPerParticle) {
        this.energyPerParticle = energyPerParticle;
    }

    public double getDenergyPerParticle() {
        return denergyPerParticle;
    }

    public void setDenergyPerParticle(double denergyPerParticle) {
        this.denergyPerParticle = denergyPerParticle;
    }

    public double getPressure() {
        return pressure;
    }

    public void setPressure(double pressure) {
        this.pressure = pressure;
    }

    public double getNumberDensity() {
        return numberDensity;
    }

    public void setNumberDensity(double numberDensity) {
        this.numberDensity = numberDensity;
    }
 
    public double getMassPotential() {
        return massPotential;
    }

    public void setMassPotential(double massPotential) {
        this.massPotential = massPotential;
    }

    public double getDpressure() {
        return dpressure;
    }

    public void setDpressure(double dpressure) {
        this.dpressure = dpressure;
    }

    public double getDnumberDensity() {
        return dnumberDensity;
    }

    public void setDnumberDensity(double dnumberDensity) {
        this.dnumberDensity = dnumberDensity;
    }

    public double getdMassPotential() {
        return dMassPotential;
    }

    public void setdMassPotential(double dMassPotential) {
        this.dMassPotential = dMassPotential;
    }

    public double getLameLambda() {
        return lameLambda;
    }

    public void setLameLambda(double lameLambda) {
        this.lameLambda = lameLambda;
    }

    public double getShearModulus() {
        return shearModulus;
    }

    public void setShearModulus(double shearModulus) {
        this.shearModulus = shearModulus;
    }
}
