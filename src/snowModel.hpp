#ifndef SNOWMODEL_H
#define SNOWMODEL_H

namespace SnowSimulator {

typedef struct SnowModel {
  double criticalCompression;
  double criticalStretch;
  double hardeningCoefficient;
  double initialDensity;
  double initialYoungsModulus;
  double poissonRatio;
  double initialLambda;
  double initialMu;

  // Suggested values from the paper

  // SnowModel()
  //     : criticalCompression(2.5e-2), criticalStretch(7.5e-3),
  //       hardeningCoefficient(10.0), initialDensity(4.0e2),
  //       initialYoungsModulus(1.4e5), poissonRatio(0.2) {
  SnowModel()
      : criticalCompression(2.5e-2), criticalStretch(7.5e-3),
        hardeningCoefficient(20.0), initialDensity(4.0e2),
        initialYoungsModulus(1.4e5), poissonRatio(0.2) {
    initialLambda = initialYoungsModulus * poissonRatio / (1 + poissonRatio) /
                    (1 - 2 * poissonRatio);
    initialMu = initialYoungsModulus / (2 * (1 + poissonRatio));
  }
} SnowModel;

} // namespace SnowSimulator

#endif // SNOWMODEL_H
