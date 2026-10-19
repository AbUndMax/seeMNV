package com.github.abundmax.seemnv.core.fieldArray;

import java.util.Collections;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

public enum FieldAggregationStrategy {

    ELEM_MIN("elem_min") {
        @Override
        public <T extends Comparable<? super T>> T[] aggregate(T[][] values) {
            T[] result = values[0].clone();

            for (int col = 0; col < result.length; col++) {
                for (int row = 1; row < values.length; row++) {
                    if (values[row][col].compareTo(result[col]) < 0) {
                        result[col] = values[row][col];
                    }
                }
            }

            return result;
        }
    },

    ELEM_MAX("elem_max") {
        @Override
        public <T extends Comparable<? super T>> T[] aggregate(T[][] values) {
            T[] result = values[0].clone();

            for (int col = 0; col < result.length; col++) {
                for (int row = 1; row < values.length; row++) {
                    if (values[row][col].compareTo(result[col]) > 0) {
                        result[col] = values[row][col];
                    }
                }
            }

            return result;
        }
    };

    private static final Map<String, FieldAggregationStrategy> BY_LABEL;

    static {
        Map<String, FieldAggregationStrategy> strategies = new HashMap<>();

        for (FieldAggregationStrategy strategy : values()) {
            strategies.put(strategy.label, strategy);
        }

        BY_LABEL = Collections.unmodifiableMap(strategies);
    }

    private final String label;

    FieldAggregationStrategy(String label) {
        this.label = label;
    }


    /**
     * Input is an Array of values for a single sample over multiple Entries aka. FieldArrays that should be aggregated
     *
     * Sample1entry1 = [[5, 3]
     * Sample1entry2 =  [2, 5]
     * Sample1entry3 =  [4, 4]]
     * ---------------------------
     * Sample1result =  [2, 3]
     *
     * @param values
     * @return
     * @param <T>
     */
    public abstract <T extends Comparable<? super T>> T[] aggregate(T[][] values);

    public static FieldAggregationStrategy fromString(String label) {
        if (label == null) {
            throw new IllegalArgumentException(
                    "Aggregation strategy must not be null"
            );
        }

        FieldAggregationStrategy strategy =
                BY_LABEL.get(label.toLowerCase(Locale.ROOT));

        if (strategy == null) {
            throw new IllegalArgumentException(
                    "Unknown aggregation strategy: " + label
            );
        }

        return strategy;
    }
}
