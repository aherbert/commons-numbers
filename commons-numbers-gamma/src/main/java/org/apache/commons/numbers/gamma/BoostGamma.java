/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//  Copyright John Maddock 2006-7, 2013-20.
//  Copyright Paul A. Bristow 2007, 2013-14.
//  Copyright Nikhar Agrawal 2013-14
//  Copyright Christopher Kormanyos 2013-14, 2020
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

package org.apache.commons.numbers.gamma;

/**
 * Implementation of the
 * <a href="http://mathworld.wolfram.com/RegularizedGammaFunction.html">Regularized Gamma functions</a> and
 * <a href="https://mathworld.wolfram.com/IncompleteGammaFunction.html">Incomplete Gamma functions</a>.
 *
 * <p>This code has been adapted from the <a href="https://www.boost.org/">Boost</a>
 * {@code c++} implementation {@code <boost/math/special_functions/gamma.hpp>}.
 * All work is copyright to the original authors and subject to the Boost Software License.
 *
 * @see
 * <a href="https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/sf_gamma.html">
 * Boost C++ Gamma functions</a>
 */
final class BoostGamma {
    /** Private constructor. */
    private BoostGamma() {
        // intentionally empty.
    }

    // Code ported from Boost 1.77.0
    //
    // boost/math/special_functions/gamma.hpp
    // boost/math/special_functions/detail/igamma_large.hpp
    //
    // Original code comments, including measured deviations, are preserved.
    //
    // Changes to the Boost implementation:
    // - Update method names to replace underscores with camel case
    // - Explicitly inline the polynomial function evaluation
    //   using Horner's method (https://en.wikipedia.org/wiki/Horner%27s_method)
    //
    // Note:
    // Constants using the 'f' suffix are machine
    // representable as a float, e.g.
    // assert 0.0891314744949340820313f == 0.0891314744949340820313;
    // The values are unchanged from the Boost reference.
}
