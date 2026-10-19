package com.github.abundmax.seemnv.core.fieldArray;

import java.util.Arrays;
import java.util.function.IntFunction;

/**
 * A FieldArray represents all values of all samples of ONE VCF entry (line) for ONE genotypeField (e.g. GT, DP, etc..)
 *
 *
 * @param <T>
 */
public class FieldArray<T extends Comparable<? super T>> {

    public final FieldConfig fieldConfig;
    public final T[][] values;
    public final int length;

    public final IntFunction<T[]> rowFactory;
    public final IntFunction<T[][]> matrixFactory;


    public FieldArray(FieldConfig fieldConfig, IntFunction<T[]> rowFactory, IntFunction<T[][]> matrixFactory, T[][] values){
        this.fieldConfig = fieldConfig;
        this.values = values;
        this.length = values.length;
        this.rowFactory = rowFactory;
        this.matrixFactory = matrixFactory;
    }

    public FieldArray(FieldConfig fieldConfig, IntFunction<T[]> rowFactory, IntFunction<T[][]> matrixFactory, String[] initialStringArray){
        this.fieldConfig = fieldConfig;
        this.rowFactory = rowFactory;
        this.matrixFactory = matrixFactory;
        this.values = parse(initialStringArray);
        this.length = values.length;
    }

    private T[][] parse(String[] initialArray) {
        T[][] parsedArray = matrixFactory.apply(initialArray.length);

        for (int i = 0; i < initialArray.length; i++) {

            // If there is more than 1 value given at the Field, it has to be split
            String[] stringFieldArray = initialArray[i].split(",");

            parsedArray[i] = rowFactory.apply(stringFieldArray.length);

            for (int j = 0; j < stringFieldArray.length; j++) {
                parsedArray[i][j] = fieldConfig.FIELD_TYPE.parseAndCast(stringFieldArray[j]);
            }
        }

        return parsedArray;
    }

    /**
     * Combines several FieldArrays of the same type using the specified combiner
     * @param fieldArrays
     * @return
     */
    @SafeVarargs
    public static <T extends Comparable<? super T>> FieldArray<T> combine(
            FieldArray<T>... fieldArrays
    ) {
        if (fieldArrays == null || fieldArrays.length == 0) {
            throw new IllegalArgumentException("At least one FieldArray is required");
        }

        FieldArray<T> first = fieldArrays[0];
        FieldConfig fieldConfig = first.fieldConfig;
        FieldAggregationStrategy strategy = fieldConfig.AGGREGATION_STRATEGY;
        IntFunction<T[][]> matrixFactory = first.matrixFactory;

        validateCompatible(fieldArrays);

        int length = first.length;
        T[][] combinedArray = Arrays.copyOf(first.values, first.length);

        for (int sample = 0; sample < length; sample++) {
            T[][] valuesOverArraysAtSample = matrixFactory.apply(fieldArrays.length);

            for (int fieldArrayIndex = 0; fieldArrayIndex < fieldArrays.length; fieldArrayIndex++) {
                valuesOverArraysAtSample[fieldArrayIndex] =
                        fieldArrays[fieldArrayIndex].values[sample];
            }

            combinedArray[sample] = strategy.aggregate(valuesOverArraysAtSample);
        }

        return new FieldArray<>(
                fieldConfig,
                first.rowFactory,
                first.matrixFactory,
                combinedArray
        );
    }

    private static <T extends Comparable<? super T>> void validateCompatible(
            FieldArray<T>[] fieldArrays
    ) {
        int length = fieldArrays[0].length;
        FieldConfig config = fieldArrays[0].fieldConfig;

        for (FieldArray<T> fieldArray : fieldArrays) {
            if (fieldArray.length != length) {
                throw new IllegalStateException("FieldArrays have to be the same size");
            }

            if (fieldArray.fieldConfig != config) {
                throw new IllegalArgumentException("FieldArrays must have the same FieldConfig");
            }
        }
    }
}
