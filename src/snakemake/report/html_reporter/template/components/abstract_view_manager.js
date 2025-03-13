class AbstractViewManager {
    handleSelectedResult(resultPath) {
        let entry = results[resultPath];
        const mimeType = getResultMimeType(resultPath);

        switch (mimeType) {
            case "image/svg+xml":
            case "image/png":
            case "image/jpeg":
                return this.handleImg(entry, resultPath);
            case "text/html":
                return this.handleHtml(entry, resultPath);
            case "application/pdf":
                return this.handlePdf(entry, resultPath);
            default:
                return this.handleDefault(entry, resultPath);
        }
    }

    handleImg(entry, resultPath) {
        throw new Error("Unimplemented!");
    }

    handleHtml(entry, resultPath) {
        throw new Error("Unimplemented!");
    }

    handlePdf(entry, resultPath) {
        throw new Error("Unimplemented!");
    }

    handleDefault(entry, resultPath) {
        throw new Error("Unimplemented!");
    }
}