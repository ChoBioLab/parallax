class ChannelUtils {

    static def validateAndJoin(channel1, channel2, joinKey = 0) {
        return channel1
            .map { tuple ->
                if (tuple.size() < 2) {
                    error "Invalid channel structure: expected [meta, ...], got ${tuple}"
                }
                return tuple
            }
            .join(
                channel2.map { tuple ->
                    if (tuple.size() < 2) {
                        error "Invalid channel structure: expected [meta, ...], got ${tuple}"
                    }
                    return tuple
                },
                by: joinKey,
                failOnMismatch: true
            )
    }

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
}
