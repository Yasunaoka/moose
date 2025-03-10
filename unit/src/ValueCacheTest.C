//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ValueCache.h"

#include <gtest/gtest.h>

TEST(ValueCacheTest, empty)
{
  ValueCache<int> cache(2);

  int out;
  Real d2;
  EXPECT_EQ(cache.guess({1, 1}, out, d2), false);
}

TEST(ValueCacheTest, guess)
{
  ValueCache<std::vector<Real>> cache(3);
  cache.insert({4.2, 2.7, 1.6}, {1, 2});
  cache.insert({3.1, 1.2, 0.1}, {3, 4});
  cache.insert({5.8, 1.9, 3.6}, {5, 6});

  std::vector<Real> out;
  Real d2;

  EXPECT_EQ(cache.guess({0, 0, 0}, out, d2), true);
  EXPECT_EQ(out[0], 3);
  EXPECT_EQ(out[1], 4);
  EXPECT_NEAR(d2, 3.1 * 3.1 + 1.2 * 1.2 + 0.1 * 0.1, 1e-9);
}

TEST(ValueCacheTest, kNN)
{
  ValueCache<std::vector<Real>> cache(3);

  cache.insert({4.2, 2.7, 1.6}, {1, 2});
  cache.insert({3.1, 1.2, 0.1}, {3, 4});

  auto nearest_neighbors = cache.kNearestNeighbors({0, 0, 0}, 3);
  EXPECT_EQ(nearest_neighbors.size(), 2);

  cache.insert({5.8, 1.9, 3.6}, {5, 6});
  cache.insert({7.1, 2.2, 4.1}, {7, 8});
  cache.insert({2.8, 2.1, 0.6}, {9, 10});

  nearest_neighbors = cache.kNearestNeighbors({0, 0, 0}, 3);

  // First nearest neighbor
  auto [key, value, distance] = nearest_neighbors[0];
  EXPECT_NEAR(key[0], 3.1, 1e-9);
  EXPECT_NEAR(key[1], 1.2, 1e-9);
  EXPECT_NEAR(key[2], 0.1, 1e-9);
  EXPECT_EQ(value[0], 3);
  EXPECT_EQ(value[1], 4);
  EXPECT_NEAR(distance, 3.1 * 3.1 + 1.2 * 1.2 + 0.1 * 0.1, 1e-9);

  // Second nearest neighbor
  std::tie(key, value, distance) = nearest_neighbors[1];
  EXPECT_NEAR(key[0], 2.8, 1e-9);
  EXPECT_NEAR(key[1], 2.1, 1e-9);
  EXPECT_NEAR(key[2], 0.6, 1e-9);
  EXPECT_EQ(value[0], 9);
  EXPECT_EQ(value[1], 10);
  EXPECT_NEAR(distance, 2.8 * 2.8 + 2.1 * 2.1 + 0.6 * 0.6, 1e-9);
}

TEST(ValueCacheTest, rebuildTree)
{
  ValueCache<int> cache(1);
  for (int i = 0; i < 110; ++i)
    cache.insert({Real(i)}, i);

  int out;
  Real d2;
  EXPECT_EQ(cache.guess({57.3}, out, d2), true);
  EXPECT_EQ(out, 57);
  EXPECT_NEAR(d2, 0.3 * 0.3, 1e-9);
}
