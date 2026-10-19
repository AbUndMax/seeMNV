package com.github.abundmax.seemnv.core.fieldArray;

import java.util.function.Function;
import java.util.function.IntFunction;

public class ArrayCombinator <T extends Comparable<? super T>> {


    private FieldAggregationStrategy strategy;
    private final IntFunction<T[][]> matrixFactory;

    public ArrayCombinator(Function<T[][], T[]> combinationFunction,
                           IntFunction<T[][]> matrixFactory) {
        this.combinationFunction = combinationFunction;
        this.matrixFactory = matrixFactory;
    }


    @SafeVarargs
    public final T[][] combine(FieldArray<T>... fieldArrays) {

        int length = fieldArrays[0].length;

        // FieldArrays of one type are assumed equal length! (otherwise a combination is not possible!)
        for (int i = 0; i < fieldArrays.length; i++) {
            if (fieldArrays[i].length != length) {
                throw new IllegalStateException("FieldArrays has to be the same size");
            }
        }

        // [[Sample1values], [Sample2Vvalues], [sample3values]...]
        T[][] combinedArray = matrixFactory.apply(length);

        // For each sample, loop over all arrays and collect their value for that sample. Then combine
        // the values to one "common" value that suffices the aggregation criteria
        for (int sample = 0; sample < length; sample++) {

            T[][] valuesOverArraysAtSample = matrixFactory.apply(fieldArrays.length);

            for (int fieldArray = 0; fieldArray < fieldArrays.length; fieldArray++) {
                valuesOverArraysAtSample[fieldArray] = fieldArrays[fieldArray].values[sample];
            }

            combinedArray[sample] = strategy.aggregate(valuesOverArraysAtSample);
        }

        return combinedArray;
    }


}
