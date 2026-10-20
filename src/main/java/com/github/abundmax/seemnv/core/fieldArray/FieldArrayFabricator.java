package com.github.abundmax.seemnv.core.fieldArray;

import java.util.Map;
import java.util.function.Function;

import static com.github.abundmax.seemnv.core.fieldArray.FieldType.*;


public class FieldArrayFabricator {


    private static final Map<FieldType, Function<FieldConfig, Function<String[], FieldArray<?>>>> REGISTRY = Map.of(
            INTEGER, fieldConfig -> initialValues -> new FieldArray<Integer>(
                    fieldConfig,
                    Integer[]::new,
                    Integer[][]::new,
                    initialValues
            ),
            FLOAT, fieldConfig -> initialValues -> new FieldArray<>(
                    fieldConfig,
                    Float[]::new,
                    Float[][]::new,
                    initialValues
            ),
            BOOLEAN, fieldConfig -> initialValues -> new FieldArray<>(
                    fieldConfig,
                    Boolean[]::new,
                    Boolean[][]::new,
                    initialValues
            ),
            STRING, fieldConfig -> initialValues -> new FieldArray<>(
                    fieldConfig,
                    String[]::new,
                    String[][]::new,
                    initialValues
            )
    );

    public static Function<String[], FieldArray<?>> get(FieldConfig fieldConfig) {
        Function<FieldConfig, Function<String[], FieldArray<?>>> factory = REGISTRY.get(fieldConfig.fieldType);

        if (factory == null) {
            throw new IllegalArgumentException(
                    "Unsupported field type: " + fieldConfig.fieldType.getLabel()
            );
        }

        return factory.apply(fieldConfig);
    }
}
