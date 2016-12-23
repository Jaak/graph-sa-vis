#pragma once

// TODO: this should all be dynamic, ofc

// Distances given in centimeters:
constexpr float BoundingBoxWidth = 30;
constexpr float BoundingBoxHeight = 20;
constexpr float BoxPadding = 0.5;

// SA Weights:
constexpr float CrossingWeight = 0.3;
constexpr float ClosenessWeight = 100;
constexpr float EdgeDistanceWeight = 1;

constexpr float EdgeFromVertexDistanceWeight = 0;
constexpr float AngleWeight = 0;

// SA parameters:
constexpr float GammaFactor = 0.98;
constexpr int StepFactor = 200;
constexpr float InitialTemperature = 1e6;
constexpr float MinimumMoveRadius = BoundingBoxHeight / 200;