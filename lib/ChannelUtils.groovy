class ChannelUtils {

    /**
     * Simple utility for standardizing metadata across channels
     * This is the only utility that adds real value without complexity
     */
    static def standardizeMeta(channel) {
        return channel.map { meta, files ->
            def std_meta = [
                id: meta.id ?: error("Missing required field: id"),
                sample: meta.sample ?: meta.id,
                condition: meta.condition ?: 'control',
                batch: meta.batch ?: 'batch1'
            ]
            [std_meta, files]
        }
    }

    /**
     * Basic validation helper - only for critical checks
     */
    static def validateMeta(meta) {
        if (!meta || !meta.id) {
            error "Invalid metadata: missing required 'id' field"
        }
        return meta
    }
}
