package com.github.abundmax.seemnv.core.fieldArray;

import java.util.function.Function;
import java.util.function.IntFunction;

public class FieldValueParser<T> {

    private final Function<String, T> parserFunction;
    private final IntFunction<T[]> rowFactory;
    private final IntFunction<T[][]> matrixFactory;

    public FieldValueParser(Function<String, T> parserFunction,
                            IntFunction<T[]> rowFactory,
                            IntFunction<T[][]> matrixFactory) {
        this.parserFunction = parserFunction;
        this.rowFactory = rowFactory;
        this.matrixFactory = matrixFactory;
    }

    public T[][] parse(String[] initialArray) {
        T[][] parsedArray = matrixFactory.apply(initialArray.length);

        for (int i = 0; i < initialArray.length; i++) {

            // If there is more than 1 value given at the Field, it has to be split
            String[] fieldArray = initialArray[i].split(",");

            parsedArray[i] = rowFactory.apply(fieldArray.length);

            for (int j = 0; j < fieldArray.length; j++) {
                parsedArray[i][j] = parserFunction.apply(fieldArray[j].trim());
            }
        }

        return parsedArray;
    }
}
