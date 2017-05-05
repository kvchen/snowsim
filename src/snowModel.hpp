struct SnowModel {
  double criticalCompression;
  double criticalStretch;
  double hardeningCoefficient;
  double initialDensity;
  double initialYoungsModulus;
  double poissonRatio;

  // Suggested values from the paper

  SnowModel()
      : criticalCompression(2.5e-2), criticalStretch(7.5e-3),
        hardeningCoefficient(10.0), initialDensity(4.0e2),
        initialYoungsModulus(1.4e5), poissonRatio(0.2) {}
};
