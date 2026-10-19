package com.github.abundmax.seemnv.core.fieldArray;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public enum FieldType {

    STRING("string", String.class) {
        @Override
        public String parseValue(String value) {
            return String.valueOf(value);
        }
    },

    INTEGER("integer", Integer.class) {
        @Override
        public Integer parseValue(String value) {
            return Integer.valueOf(value);
        }
    },

    FLOAT("float", Float.class) {
        @Override
        public Float parseValue(String value) {
            return Float.valueOf(value);
        }
    },

    BOOLEAN("boolean", Boolean.class) {
        @Override
        public Boolean parseValue(String value) {
            return Boolean.valueOf(value);
        }
    };

    private static final Map<String, FieldType> BY_LABEL;

    static {
        Map<String, FieldType> fieldTypes = new HashMap<>();

        for (FieldType type : values()) {
            fieldTypes.put(type.label, type);
        }

        BY_LABEL = Collections.unmodifiableMap(fieldTypes);
    }

    private final String label;
    private final Class<?> type;

    FieldType(String label, Class<?> type) {
        this.label = label;
        this.type = type;
    }

    public abstract Object parseValue(String value);

    @SuppressWarnings("unchecked")
    public <T> T parseAndCast(String value) {
        return (T) parseValue(value.trim());
    }

    public String getLabel() {
        return label;
    }

    public Class<?> getType() {
        return type;
    }

    public static FieldType fromString(String label) {
        FieldType fieldType = BY_LABEL.get(label.toLowerCase());

        if (fieldType == null) {
            throw new IllegalArgumentException("Unsupported field type: " + label);
        }

        return fieldType;
    }
}
