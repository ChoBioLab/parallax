class ChannelUtils {

    static def validateAndJoin(channel1, channel2, joinKey = 0, strict = true) {
        // Validate input channels first
        def validated_ch1 = channel1
            .map { tuple ->
                validateTuple(tuple, "channel1")
                return tuple
            }

        def validated_ch2 = channel2
            .map { tuple ->
                validateTuple(tuple, "channel2")
                return tuple
            }

        // Perform join with error handling
        return validated_ch1
            .join(
                validated_ch2,
                by: joinKey,
                failOnMismatch: strict,
                remainder: !strict
            )
            .map { tuple ->
                // Post-join validation
                if (tuple.size() < 3) {
                    error "Join failed: expected at least 3 elements after join, got ${tuple.size()}"
                }
                return tuple
            }
    }

    static def safeJoin(channel1, channel2, joinKey = 0) {
        // Create indexed channels for tracking
        def indexed_ch1 = channel1
            .map { tuple -> 
                def meta = tuple[0]
                def key = getJoinKey(meta, joinKey)
                [key, tuple]
            }

        def indexed_ch2 = channel2
            .map { tuple ->
                def meta = tuple[0]
                def key = getJoinKey(meta, joinKey)
                [key, tuple]
            }

        // Join and validate
        return indexed_ch1
            .join(indexed_ch2, failOnMismatch: true)
            .map { key, tuple1, tuple2 ->
                // Merge tuples, keeping first meta
                def result = [tuple1[0]] + tuple1[1..-1] + tuple2[1..-1]
                validateJoinResult(result, key)
                return result
            }
    }

    static def synchronizedJoin(channels, joinKey = 0) {
        if (channels.size() < 2) {
            error "synchronizedJoin requires at least 2 channels"
        }

        def result = channels[0]
        for (int i = 1; i < channels.size(); i++) {
            result = safeJoin(result, channels[i], joinKey)
        }
        return result
    }

    static def validateTuple(tuple, channelName = "unknown") {
        if (!tuple || tuple.size() < 2) {
            error "Invalid tuple in ${channelName}: expected [meta, ...], got ${tuple}"
        }

        def meta = tuple[0]
        if (!meta || !meta.id) {
            error "Invalid meta in ${channelName}: missing required 'id' field"
        }

        return tuple
    }

    static def getJoinKey(meta, keySpec) {
        if (keySpec instanceof Integer) {
            return meta.id  // Default to ID
        } else if (keySpec instanceof String) {
            return meta[keySpec] ?: error("Missing join key '${keySpec}' in meta: ${meta}")
        } else {
            return keySpec.call(meta)
        }
    }

    static def validateJoinResult(result, key) {
        if (!result || result.size() < 2) {
            error "Join validation failed for key '${key}': invalid result ${result}"
        }
        return result
    }

    static def standardizeMeta(channel) {
        return channel.map { meta, files ->
            def std_meta = [
                id: meta.id ?: error("Missing required field: id"),
                sample: meta.sample ?: meta.id,
                condition: meta.condition ?: 'control',
                batch: meta.batch ?: 'batch1',
                timestamp: System.currentTimeMillis()
            ]
            [std_meta, files]
        }
    }

    static def withTimeout(channel, timeoutSeconds = 300) {
        return channel
            .map { tuple ->
                def startTime = System.currentTimeMillis()
                // Add timeout tracking to meta
                def meta = tuple[0]
                meta.processing_start = startTime
                return tuple
            }
    }

    static def collectWithValidation(channel, expectedCount = null) {
        def collected = channel.collect()

        if (expectedCount != null) {
            return collected.map { items ->
                if (items.size() != expectedCount) {
                    error "Collection validation failed: expected ${expectedCount} items, got ${items.size()}"
                }
                return items
            }
        }

        return collected
    }
}
